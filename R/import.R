

#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @noRd
NULL





#' @importFrom dplyr arrange bind_cols bind_rows case_when count distinct everything filter group_by lag lead left_join mutate mutate_all n near pull rename select slice summarise tally ungroup
#' @importFrom ggplot2 aes facet_wrap geom_line geom_rect ggplot scale_y_continuous theme_bw xlab
#' @importFrom glue glue
#' @importFrom purrr flatten keep map map_dbl map_dfr pmap
#' @importFrom readr write_tsv
#' @importFrom rlang .data
#' @importFrom stats rpois runif setNames
#' @importFrom stringr str_c
#' @importFrom tibble as_tibble enframe is_tibble tibble
#' @importFrom tidyr nest pivot_longer pivot_wider separate unite unnest
#' @importFrom utils write.table
NULL



# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "1",
      "2",
      "bp",
      "BY",
      "GSP_opts",
      "abs_column",
      "alle_int",
      "allele",
      "chr_xmax",
      "chr_xmin",
      "chr_ymax",
      "chr_ymin",
      "chrom",
      "chrom_f",
      "chrom_len",
      "data",
      "end",
      "end_pos",
      "gam_tibs",
      "gamete_index",
      "gamete_segments",
      "gene_copy",
      "gpp",
      "group",
      "group_length",
      "group_origin",
      "gs_column",
      "gsp_init",
      "h",
      "hap1",
      "hap2",
      "haplo",
      "id",
      "idx",
      "ind",
      "index",
      "indiv",
      "link_pos",
      "locus",
      "m_subscript_matrix",
      "ma",
      "map_stuff",
      "matrix_row",
      "max_end_pos",
      "max_rec",
      "mid_pos",
      "morgans",
      "names_col",
      "new_tib",
      "next_start",
      "pa",
      "par1",
      "par2",
      "ped_sample_id",
      "ped_samples",
      "pheno",
      "pop",
      "pop_mat",
      "pop_origin",
      "pos",
      "rec_prob",
      "rows",
      "rs_founder_haplo",
      "s_col",
      "samp_index",
      "segged",
      "sex_code",
      "start",
      "start_pos",
      "tmp_seg_names",
      "tot_length",
      "unit",
      "useit",
      "variant_id",
      "yval",
      "z_col"
    )
  )
}
