% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gsp_for_validity_and_saturation.R
\name{check_gsp_for_validity_and_saturation}
\alias{check_gsp_for_validity_and_saturation}
\title{Check that GSP does not reuse genetic material and yet uses all of it}
\usage{
check_gsp_for_validity_and_saturation(GP)
}
\arguments{
\item{GP}{A gsp in list format as produced by the function
\code{prep_gsp_for_hap_dropping()}.  See the documentation for the return
object of \code{prep_gsp_for_hap_dropping()} for a description.}
}
\value{
This function does not return anything.
}
\description{
These conditions can be checked for a GSP with no inbreeding loops
simply by ensuring that the amount of genetic material coming into
each individual is the same as the amount going out (either as
segregated gametes or as samples).  If the amount of material coming
out of any individual in the GSP is greater than the amount coming
in, then an error is thrown.  If the amount coming out is less than
the amount coming in, then a warning about the GSP is thrown.
Messages printed via \code{message()} and \code{warning()} indicate which
individuals in the GSP are problematic.  All problematic individuals
are listed before an error is thrown with \code{stop()}.
}
\examples{
# get the 13 member pedigree in tibble form as the package
# and turn it into a list
GP <- prep_gsp_for_hap_dropping(GSP)

# check it. (This passes)
check_gsp_for_validity_and_saturation(GP)

## The following will show a failure, so we wrap it in tryCatch
## so CRAN check does not flag it as a problem.
# Read in a gsp with errors and then make sure all the
# error in it are caught
bad <- readr::read_csv(system.file("extdata/13-member-ped-with-errors.csv", package = "gscramble"))

# check_gsp_for_validity_and_saturation() is called internally from
# within prep_gsp_for_hap_dropping(), after creating a list-from GSP.
# This will show the error produced by check_gsp_for_validity_and_saturation().
badL <- tryCatch(
  prep_gsp_for_hap_dropping(bad),
  error = function(x) 0,
  warning = function(x) 0
)


}
