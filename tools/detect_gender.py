from json import dump
from os.path import basename, dirname, join, realpath
from pprint import pprint

from pandas import read_csv
from numpy import nanmedian, nanmean, arange, NaN
from pysam import TabixFile, asTuple

GENOME_APP_DIRECTORY_PATH = dirname(dirname(realpath(__file__)))

GENOME_APP_NAME = basename(GENOME_APP_DIRECTORY_PATH)

INPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'input')
PERSON_DIRECTORY_PATH = join(INPUT_DIRECTORY_PATH, 'person')
GRCH_DIRECTORY_PATH = join(INPUT_DIRECTORY_PATH, 'grch')

TOOLS_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'tools')
OUTPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'output')
MEDIA_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'media')

REGION_FILE = join(INPUT_DIRECTORY_PATH, 'non_PAR_region.bed')

VCF_FILE = join(PERSON_DIRECTORY_PATH, 'genome.vcf.gz')

tbx = TabixFile(VCF_FILE)


def get_format_index():
    """
    Get the index of the FORMAT field in VCF file
    Returns:
        int: index of FORMAT field
    """
    try:
        for row in tbx.header:
            line = row.decode('UTF-8')
            if line.startswith('#CHROM') and 'FORMAT' in line:
                index = line.split('\t').index('FORMAT')
                return index
    except NameError:
        print('FORMAT field not found')


def get_format_field_index(field):
    """
    Get the index of the specified field in the FORMAT string
    Arguments:
        field (str): field to extract
    Returns:
        int: index of field in FORMAT string
    """
    try:
        for row in tbx.fetch():
            field_index = row.split('\t')[get_format_index()].split(':').index(
                field)
            return field_index
    except ValueError:
        print('{} not found'.format(field))


def get_field_value(row, format_index, field_index):
    """
    Get the index of the specified field in the FORMAT string
    Arguments:
        row (str): field to extract
        format_index (int): index of FORMAT string in VCF files
        field_index (int): field to field to extract in FORMAT string
    Returns:
        str: value of specified field
    """
    try:
        return row[format_index + 1].split(':')[field_index]
    except IndexError:
        print('Invalid index')


def check_chr_format():
    """
    Check the format of the CHROM field in VCF file
    Returns:
        bool: boolean of whether format is 'chr1' or '1'
    """
    try:
        for row in tbx.header:
            line = row.decode('UTF-8')
            if line.startswith('##contig'):
                return 'chr' in line.split('ID=')[1]
    except ValueError:
        print('Contig in header not in correct format.')


def get_info_from_variants(chrom, start, stop, field):
    """
    Check the format of the CHROM field in VCF file
    Arguments:
        chrom (str): chromosome of search region
        start (str): start position of search region
        stop (str): stop position of search region
        field (str): field to extract from VCF file
    Returns:
        list (tuple): list of tuples by (chromosome, field value)
    """
    tuples = []
    format_index = get_format_index()
    field_index = get_format_field_index(field)
    try:
        for entry in tbx.fetch(
                chrom, int(start), int(stop), parser=asTuple()):

            if float(entry[5]) > 20:
                tuples.append((chrom, get_field_value(entry, format_index,
                                                      field_index)))
    except ValueError:
        print(
            "No variants found in region {}:{}:{}".format(chrom, start, stop))
        pass

    return tuples


def get_variant_info(region, field):
    """
    Build the dictionary of specified field for all regions specified in BED file
    Arguments:
        field (str): field to extract from VCF file
    Returns:
        dict: {'chromosome': [field value1, field value2, ...]}
    """
    variant_info = {}

    for chrom, start, stop in zip(region.CHROM, region.START, region.STOP):

        # Check chromosome format
        if check_chr_format():
            chrom = 'chr' + chrom

        # Get variant info per region
        variants = get_info_from_variants(chrom, start, stop, field)

        # Build variant dict
        for entry in variants:
            try:
                if entry[0] not in variant_info:
                    variant_info[entry[0]] = [entry[1]]
                else:
                    variant_info[entry[0]].append(entry[1])
            except ValueError:
                pass

    return variant_info


def detect_gender_from_read_depth(region, threshold=0.2):
    """
    Detect gender by comparing the read depth of the sex chromosomes against that of autosomal chromosomes
    Arguments:
        threshold (double): threshold difference between read depths
    Returns:
        str: gender
    """
    field = 'DP'

    # Get field for each variant
    variant_info = get_variant_info(region, field)

    # Convert read depth string to int
    for key in variant_info:
        variant_info[key] = [
            int(x) if x.isdigit() else NaN for x in variant_info[key]
        ]

    # Get correct chromosome list format
    if check_chr_format():
        chromosomes = ['chr' + str(x) for x in list(arange(1, 23))]
        chrx = 'chrX'
        chry = 'chrY'
    else:
        chromosomes = [str(x) for x in list(arange(1, 23))]
        chrx = 'X'
        chry = 'Y'

    # Calculate average coverage and coverage over sex chromosomes
    average_coverage = nanmean(
        [nanmedian(variant_info[x]) for x in variant_info if x in chromosomes])
    x_coverage = nanmedian(variant_info[chrx])

    try:
        y_coverage = nanmedian(variant_info[chry])
    except KeyError:
        y_coverage = 0
        print('No coverage found for Y chromosome.')

    # Determine gender based on coverage in X and Y chromosomes
    if (x_coverage < average_coverage * (1 - threshold) and
            x_coverage > average_coverage * threshold and
            y_coverage < average_coverage * (1 - threshold)):
        gender = 'male'
    elif (y_coverage < average_coverage * threshold and
          x_coverage > average_coverage * (1 - threshold) and
          x_coverage < average_coverage * (1 + threshold)):
        gender = 'female'
    else:
        gender = None

    return gender


def detect_gender_from_het_hom_ratio(region, threshold=0.2):
    """
    Detect gender by comparing the het/hom ratio of the sex chromosomes against that of autosomal chromosomes
    Arguments:
        threshold (double): threshold difference between het/hom ratios
    Returns:
        str: gender
    """
    field = 'GT'
    het_hom_ratio = {}

    # Get field for each variant
    variant_info = get_variant_info(region, field)

    for key in sorted(variant_info.keys()):
        hom = len([x for x in variant_info[key] if x == '1/1' or x == '2/2'])
        het = len([
            x for x in variant_info[key]
            if x == '0/1' or x == '1/0' or x == '1/2' or x == '2/1' or
            x == '0/2' or x == '2/0'
        ])

        het_hom_ratio[key] = het / hom

    # Get correct chromosome list format
    if check_chr_format():
        chromosomes = ['chr' + str(x) for x in list(arange(1, 23))]
        chrx = 'chrX'
    else:
        chromosomes = [str(x) for x in list(arange(1, 23))]
        chrx = 'X'

    average_ratio = nanmean(
        [het_hom_ratio[x] for x in het_hom_ratio if x in chromosomes])
    x_ratio = het_hom_ratio[chrx]

    if x_ratio < average_ratio * (1 - threshold):
        gender = 'male'
    elif x_ratio > average_ratio * (1 - threshold):
        gender = 'female'

    return gender


def detect_gender():
    """
    Detect gender using 2 methods
    1) Using read depth
    2) Using het/hom ratio
    Will always return a gender
    Returns:
        str: gender
    """

    # Remove PAR regions
    region = read_csv(REGION_FILE, sep='\t')

    print('Calculating gender from read depth...')
    gender = detect_gender_from_read_depth(region)

    if not gender:
        print(
            'Could not calculate gender from read depth. Calculating gender from het/hom ratio...'
        )
        gender = detect_gender_from_het_hom_ratio(region)

    output = {}
    output['Gender'] = gender

    # Write results as ../output/output.json
    output_json_file_path = join(OUTPUT_DIRECTORY_PATH, 'output.json')
    with open(output_json_file_path, 'w') as f:
        dump(output, f, indent=2, sort_keys=True)

    # Summarize
    print('This Genome App ran and produced {}.'.format(output_json_file_path))
    pprint(output)

    return gender
