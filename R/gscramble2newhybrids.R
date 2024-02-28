

#' Convert 'gscramble' output to newhybrids format
#'
#' This function turns character-based alleles into integers
#' and writes the necessary headers, etc.  It preferentially uses
#' the "id" column if it exists in `M$ret_ids`. Otherwise it uses
#' the `indiv` column for the sample names.
#'
#' It allows you to set the -s and -z through some regular expression mapping.
#' @param M the output from `segments2markers()` from 'gscramble'. This could
#' have an added `id` column on it, which will then be used for the
#' sample names.
#' @param M_meta the Marker meta data file.
#' @param z A vector of length two. The values
#' are regular expressions that the sample names that you want to have
#' -z 0 or -z 1 should match.  For example `c("SH", "CCT")` means any
#' sample matching "SH" would get z0 and any sample matchine "CCT" would
#' get z1.
#' @param s a single regular expressions that matches individuals that
#' should be given the -s option. For example "SH|CCT"
#' @param retain a vector of loci to retain.
#' @param outfile path to the file to write the newhybrids data set to. For CRAN
#' compliance, this is, by default, a temp file. But you can change it to be
#' anything valid.
#' @details This function relies a lot on some tidyverse functions
#' for pivoting, etc.  As such, it is not intended for data sets with
#' tens of thousands of markers.  You oughtn't be using NewHybrids with
#' so many markers, anyway!
#' @return A list with three components:
#'   - `outfile`: outfile name of saved data.
#'   - `genos`: Genotypes
#'   - `allele_names`: Allele names.
#' @export
#' @examples
#' # get output from segments2markers():
#' example("segments2markers")
#' # copy that result to a new variable
#' M <- s2m_result
#'
#' # then run it
#' gscramble2newhybrids(M, M_meta)
gscramble2newhybrids <- function(
    M,
    M_meta,
    z = NULL,
    s = NULL,
    retain = NULL,
    outfile = tempfile()
  ) {
  headers <- paste(rep(M_meta$variant_id, each = 2), "____", c(1,2), sep = "")

  if(any(names(M$ret_ids) == "id")) {
    IDs <- M$ret_ids %>% mutate(useit = id)
  } else {
    IDs <- M$ret_ids %>% mutate(useit = indiv)
  }
  geno <- M$ret_geno
  colnames(geno) <- headers
  geno_tib <- bind_cols(
    IDs %>% select(useit),
    as_tibble(geno)
  )

  geno_long <- geno_tib %>%
    pivot_longer(-useit, names_to = c("locus", "gene_copy"), values_to = "allele", names_sep = "____")

  # now, turn the alleles into integers:
  tmp <- geno_long %>%
    group_by(locus) %>%
    mutate(
      alle_int = sprintf("%02d", as.integer(as.factor(allele))),
      alle_int = ifelse(is.na(allele), "00", alle_int)
    ) %>%
    ungroup()

  # at this juncture, if requested, we will filter down the loci included
  if(!is.null(retain)) {
    tmp2 <- tmp %>%
      filter(locus %in% retain)
  } else {
    tmp2 <- tmp
  }

  # store the allele names (and counts) to return with the newhybs genos
  allele_names <- tmp2 %>%
    count(locus, allele, alle_int)

  # now, pivot the thing back, after smushing the two alleles together
  genos2 <- tmp2 %>%
    select(-allele) %>%
    pivot_wider(
      names_from = c(gene_copy),
      values_from = alle_int,
    ) %>%
    mutate(alle_int = paste0(`1`, `2`)) %>%
    select(-`1`, -`2`) %>%
    pivot_wider(names_from = locus, values_from = alle_int)



  # now, we just have to add the names, etc.
  genos3 <- genos2 %>%
    mutate(names_col = paste("n", useit), .after = useit)

  if(!is.null(z)) {
    stopifnot(length(z) == 2)
    genos4 <- genos3 %>%
      mutate(
        z_col = case_when(
          str_detect(useit, z[1]) ~ " z0",
          str_detect(useit, z[2]) ~ " z01",
          TRUE ~ ""
        ),
        .after = names_col
      )
  } else{
    genos4 <- genos3 %>%
      mutate(z_col = "", .after = names_col)
  }
  if(!is.null(s)) {
    stopifnot(length(s) == 1)
    genos5 <- genos4 %>%
      mutate(
        s_col = case_when(
          str_detect(useit, s) ~ "s",
          TRUE ~ ""
        ),
        .after = z_col
      )
  } else {
    genos5 <- genos4 %>%
      mutate(
        s_col = "",
        .after = z_col
      )
  }

  genos6 <- genos5 %>%
    mutate(opt_col = paste0(names_col, z_col, s_col), .after = names_col) %>%
    select(-z_col, -s_col, -names_col)

  # get one with index numbers...
  genos7 <- genos6 %>%
    mutate(useit = 1:n())

  # Now that is done, we just need to write it out
  num_loci <- (ncol(genos6) - 2)
  locnames <- names(genos6)[-c(1, 2)]
  cat("NumIndivs ", nrow(genos6), "\n", sep = "", file = outfile)
  cat("NumLoci ", num_loci, "\n", sep = "", file = outfile, append = TRUE)
  cat("Digits 2\nFormat Lumped\n", file = outfile, append = TRUE)
  cat("\nLocusNames ", locnames, sep = " ", file = outfile, append = TRUE)
  cat("\n\n", file = outfile, append = TRUE)

  write_tsv(genos7, col_names = FALSE, file = outfile, append = TRUE)

  # and return a list
  list(
    outfile = outfile,
    genos = genos7,
    allele_names = allele_names
  )

}
