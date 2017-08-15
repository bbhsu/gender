# Gender

Gender is the range of characteristics pertaining to, and differentiating between, masculinity and femininity. Biological gender can be determined from your DNA.

## How this app works?

This genome app uses up to 2 methods to check gender:

1) Checking the average read depth of the sex chromosomes (X, Y) compared to the autosomal chromosomes (chr1-22)
2) Checking the het/hom ratio of the sex chromosomes compared to the autosomal chromosomes

For females, the average read depth of the Y chromosome should be 0 because no Y chromosome exists - but in reality this number can be greater than 0 because current bioinformatics algorithms are not perfect. For males, the average read depth of the X and Y chromosomes should be about half of that of the autosomal chromosomes because only one copy of each exist.

For females, the het/hom ratio of the X chromosome should be around the same as that of the autosomal chromosomes because two copies of the X chromosome exist. For males, the het/hom ratio of the X chromosome should be close to 0 because the majority of genotype calls will be 1/1 (depending on the variant caller) and the het/hom ratio of the Y chromosome should be close to that of autosomal chromosomes because many regions in the Y chromosome overlap with the X chromosome.
## Want to run this app on your own DNA?

All you need is a tabixed VCF.gz file.

1. Clone this repository.
```
git clone https://github.com/bbhsu/gender
```
2. Replace 'genome.vcf.gz' in the VCF_FILE variable with the path to your tabixed VCF.gz file.

[Genome App](https://www.guardiome.com/genome-apps/) powered by [Guardiome](https://www.guardiome.com/)

<div>
  <img src="media/guardiome-logo.png" width="150" height="150">
</div>
