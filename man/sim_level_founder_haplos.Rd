% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_level_founder_haplos.R
\name{sim_level_founder_haplos}
\alias{sim_level_founder_haplos}
\title{computes the simulation-level founder haplotype index for each founder haplo}
\usage{
sim_level_founder_haplos(S)
}
\arguments{
\item{S}{tibble of segments like that produced by \code{\link{segregate}}.}
}
\value{
This function returns a result that is basically the output of \code{segregate()} with
an additional column added to it:  \code{sim_level_founder_haplo}.  This is the index
of the haplotype within each \code{group_origin} that should be used. For details of the
other columns in the output tibble, see the documentation for \code{\link{segregate}}.
}
\description{
This takes the output of \code{segregate()} and deals with the multiple gpp's and reps
to come up with a unique index for each found haplotype, so that those haplotypes
can all, eventually, be accessed easily out of the genotype matrix.
Along the way, this function does some light checking to make sure that the
\code{rs_founder_haplo} values are dense within \code{gpp} and \code{index} as they should be.
}
\examples{
#### Get output from segregate to use as input ####
# We construct an example here where we will request segregation
# down a GSP with two F1s and F1B backcrosses between two hypothetical
# populations, A and B.
gsp_f1f1b <- create_GSP("A", "B", F1 = TRUE, F1B = TRUE)

# We will imagine that in our marker data there are three groups
# labelled "grp1", "grp2", and "grp3", and we want to create the F1Bs with backcrossing
# only to grp3.
reppop <- tibble::tibble(
  index = as.integer(c(1, 1, 2, 2)),
  pop = c("A", "B", "A", "B"),
  group = c("grp3", "grp1", "grp3", "grp2")
)

# combine those into a request
request <- tibble::tibble(
  gpp = list(gsp_f1f1b),
  reppop = list(reppop)
)

# now run it through segregate()
set.seed(5)  # just for reproducibility in example...
simSegs <- segregate(request, RecRates)

#### Now we can run those through sim_level_founder_haplos() ####
fh <- sim_level_founder_haplos(simSegs)
fh
}
