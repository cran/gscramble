% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{M_meta}
\alias{M_meta}
\title{Metadata for 100 molecular markers}
\format{
A tibble with three columns: chrom, pos, and variant_id.
}
\source{
Chromosome, position, and variant IDs are from the Sus scrofa 10.2 genome \url{https://www.ncbi.nlm.nih.gov/assembly/GCF_000003025.5/}.
}
\description{
This data set contains the metadata for the 100 most divergent SNP loci for three feral swine populations
sampled in the United States. To make the dataset more computationally efficient, only 3 chromosomes were
used (12,17, and 18). The metadata for the SNP loci is in a tibble with three columns:
chrom (character), pos (numeric), variant_id (character).
The column 'chrom' contains the chromosome ID where the SNP is located, the column 'pos' gives the base pair
location on the chromosome, and the column 'variant_id' gives the name of the SNP.
Each row of the metadata tibble contains the metadata for a given SNP locus.
The individual genotypes for each of these SNP loci can be found in \code{\link{Geno}}.
}
