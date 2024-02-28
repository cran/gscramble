#' Take the output of rearrange_genos and permute everyone by population
#'
#' This is done prior to assigning random genomic fragments of individuals in the
#' sample to the founders of the GSP, to be dropped to the samples.
#' @param GS the tibble that is the output from rearrange_genos
#' @param preserve_haplotypes If TRUE then the Geno data is assumed phased
#' (first allele at an individual on one haplotype and second allele on the
#' other) and those haplotypes are preserved in this permutation of
#' genomic material amongst the founders.
#' @param preserve_individuals If TRUE then whole individuals are permuted
#' around the data set and the two gene copies at each locus are randomly
#' permuted within each individual.  If `preserve_individuals = "BY_CHROM"`,
#' then the the two copies of each chromosome in an individual are permuted
#' together.  Thus a permuted individual may have two copies of one chromosome
#' from one individual, and two copies of another chromosome from a different
#' individual.  (If `preserve_haplotypes = TRUE` then
#' the gene copies are not permuted within individuals. You should only ever
#' use `preserve_haplotypes = TRUE` if you have phased data.)
#' @return Returns a list of the same format as the output of `rearrange_genos`.
#' Plus one additional component. Each component of the return list is itself
#' an unnamed list with one component (makes it easier to use `bind_rows` to
#' create a tibble of list columns from these).  The components, once unlisted are:
#' - `G`: a matrix---the original genotype data matrix
#' - `I`: the I_meta tibble
#' - `M`: the M_meta tibble
#' - `G_permed`: the genotype matrix after permutation.
#' @export
#' @examples
#' # first get the output of rearrange_genos
#' RG <- rearrange_genos(Geno, I_meta, M_meta)
#'
#' # then permute by the populations
#' PG <- perm_gs_by_pops(RG)
perm_gs_by_pops <- function(GS, preserve_haplotypes = FALSE, preserve_individuals = FALSE) {

  row_groups <- NULL


  if((preserve_individuals == FALSE && preserve_haplotypes == TRUE) || preserve_individuals == "BY_CHROM") {
    # get a row_groups list that gives the indexes of markers on different chromosomes
    Mm = GS$M[[1]]
    row_groups <- Mm %>%
      ungroup() %>%
      mutate(idx = 1:n()) %>%
      group_by(chrom) %>%
      summarise(rows = list(idx)) %>%
      pull(rows)
  }

  scrambit <- GS$I[[1]] %>%
    group_by(group) %>%
    summarise(pop_mat = list(
      mat_scramble(
        GS$G[[1]][, abs_column],
        row_groups = row_groups,
        preserve_individuals = preserve_individuals,
        preserve_haplotypes = preserve_haplotypes
      ))) %>%
    pull(pop_mat) %>%
    do.call(what = cbind, args = .)

  GS$G_permed <- list(scrambit)

  GS
}
