#' Renumber GSP members by adding a constant to each
#'
#' This function assumes that all individuals are named as
#' numerics and that their haplotypes in hap1 and hap2
#' are named Xa and Xb,respectively, and their
#' samples are named sX, where X is an integer,
#' @param G a GSP tibble
#' @param add amount to add to each label
#' @return Returns a GSP just like the input, but with the identity numbers
#' of the individuals in it incremented by `add`.
#' @export
#' @examples
#' # get an example GSP
#' G <- create_GSP(pop1 = "p1", pop2 = "p2", F1B2 = TRUE)
#'
#'
renumber_GSP <- function(G, add) {

    mutate(G,
      ind = add + ind,
      par1 = add + par1,
      par2 = add + par2,
      hap1 = ifelse(is.na(hap1), hap1, paste0(ind, "a")),
      hap2 = ifelse(is.na(hap2), hap2, paste0(ind, "b")),
      sample = ifelse(is.na(sample), sample, paste0("s", ind))
    )
}
