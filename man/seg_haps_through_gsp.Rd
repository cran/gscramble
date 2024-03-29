% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seg_haps_through_gsp.R
\name{seg_haps_through_gsp}
\alias{seg_haps_through_gsp}
\title{Segregate haplotypes through a genome simulation pedigree list}
\usage{
seg_haps_through_gsp(
  G,
  M,
  chrom_len,
  pop_idx_sep = "--\%--",
  recomb_model = "always one"
)
}
\arguments{
\item{G}{a genome simulation pedigree as a list, ready to go}

\item{M}{map information.  Should be a tibble with start_pos end_pos
and rec_prob}

\item{chrom_len}{The length of the chromosome in base pairs}

\item{pop_idx_sep}{The string to use to separate the index of the pop
and the index of the founder for naming where chromosomal segments
come from.  It is not recommended that this be changed.}

\item{recomb_model}{Not used for now, but we will want to add that
here, eventually, to specify the kind of recombination model we are using,
and then modify recomb_point to accommodate the different choices.}
}
\description{
The main input into this is the list prepared by prep_gsp_for_hap_dropping.
This returns a list of the same form, but with a few extra variables in there,
the most important of which for future use will be the Samples, which are
lists that hold the founder segments that made it into the sampled individuals.
}
\details{
This is a low-level function that the user will not typically use.
}
\keyword{internal}
