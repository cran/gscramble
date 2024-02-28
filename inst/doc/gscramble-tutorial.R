## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gscramble)
library(tidyverse)

## -----------------------------------------------------------------------------
Geno[1:4, 1:12]

## ---- eval=FALSE--------------------------------------------------------------
#  storage.mode(IntMat) <- "character"

## -----------------------------------------------------------------------------
Geno[Geno == "-1"] <- NA

## -----------------------------------------------------------------------------
head(I_meta)

## -----------------------------------------------------------------------------
head(M_meta)

## -----------------------------------------------------------------------------
head(RecRates)

## ---- fig.align='center', out.width='50%', echo = FALSE, fig.cap="**A Simple GSP.** Here we have a pedigree with  diploid founder '1', who is formed by the union of two gametes from Population A, and diploid founder '2' with two gametes from Population B. Founder '1' and founder '2' then pass 2 gametes (red 2's) for the formation of 2 individuals represented by their descendant 'box 3'. The pink hexagon 's3' then represents the samples taken from 'box 3' with the pink 2 used to indicate the number of diploid individuals sampled."----
knitr::include_graphics("../man/figures/F1-500.png")

## -----------------------------------------------------------------------------
gspF1 <- create_GSP(pop1 = "p1", pop2 = "p2", F1 = TRUE)

gspF1

## -----------------------------------------------------------------------------
Pattern = c("Pop1", "Pop2")
RepPopSimpleF1 <- tibble(
  index = rep(1:1, each = 2),
  pop = rep(c("p1", "p2"), times = 1),
  group = Pattern
)
RepPopSimpleF1

## -----------------------------------------------------------------------------
Pattern = c("Pop1", "Pop2", "Pop1",
            "Pop3", "Pop1", "Pop4")
RepPopSimpleF1_b <- tibble(
  index = rep(1:3, each = 2),
  pop = rep(c("p1", "p2"), times = 3),
  group = Pattern
)
RepPopSimpleF1_b

## -----------------------------------------------------------------------------
Pattern = c("Pop1", "Pop2", "Pop1", "Pop3", "Pop1", "Pop4", 
            "Pop2", "Pop3", "Pop2", "Pop4", "Pop3", "Pop4")
RepPopSimpleF1_c <- tibble(
  index = rep(1:6, each = 2),
  pop = rep(c("p1", "p2"), times = 6),
  group = Pattern
)
head(RepPopSimpleF1_c)


## -----------------------------------------------------------------------------
gspComplex <- create_GSP(
  pop1 = "p1", 
  pop2 = "p2", 
  F1 = TRUE, 
  F1B = TRUE, 
  F1B2 = TRUE
)

gspComplex

## ---- fig.align='center', out.width='100%', echo = FALSE----------------------
knitr::include_graphics("../man/figures/13.svg")

## -----------------------------------------------------------------------------
Pattern = c("Pop1", "Pop2", "Pop1", 
            "Pop3", "Pop1", "Pop4")
RepPopComplex1 <- tibble(
  index = rep(1:3, each = 2),
  pop = rep(c("p1", "p2"), times = 3), 
  group = Pattern
)
RepPopComplex1

## ---- fig.align='center', out.width='100%', echo = FALSE----------------------
knitr::include_graphics("../man/figures/13-member-ped.svg")

## ----eval=FALSE---------------------------------------------------------------
#  csv <- system.file("extdata/13-member-ped.csv", package = "gscramble")
#  gsp_tib <- readr::read_csv(csv)
#  paths <- gsp2dot(g = gsp_tib, path = "images/13-member-ped")
#  # now, get rid of the dot and png files
#  file.remove(paths[1:2])

## -----------------------------------------------------------------------------
GSP

## ---- eval=FALSE--------------------------------------------------------------
#  system.file("extdata/13-member-ped.csv", package = "gscramble")

## ----eval=FALSE---------------------------------------------------------------
#  system.file("extdata/gsp4.csv", package = "gscramble")

## ---- fig.align='center', out.width='70%', echo = FALSE-----------------------
knitr::include_graphics("../man/figures/gsp4-700.png")

## -----------------------------------------------------------------------------
I_meta %>%
  count(group)

## -----------------------------------------------------------------------------
RepPop1

## -----------------------------------------------------------------------------
RepPop4

## -----------------------------------------------------------------------------
Input_tibble <- tibble(
  gpp = list(gspComplex),
  reppop = list(RepPopSimpleF1)
)

# here is what that input object looks like:
Input_tibble

## -----------------------------------------------------------------------------
set.seed(15) # for reproducibility
Segments <- segregate(
  request = Input_tibble,
  RR = RecRates
)

## -----------------------------------------------------------------------------
Segments

## -----------------------------------------------------------------------------
Segments %>%
  select(gpp:pop_origin)

## -----------------------------------------------------------------------------
Segments %>%
  select(rs_founder_haplo:group_origin)

## ---- fig.align='center', out.width='100%', echo = FALSE----------------------
knitr::include_graphics("../man/figures/13.svg")

## ---- fig.width=7.5, fig.height=9---------------------------------------------
g <- plot_simulated_chromomsome_segments(Segments, RecRates)
g

## -----------------------------------------------------------------------------
set.seed(15) # for reproducibility
# re-run segment segregation, explicitly passing it the
# marker meta data
Segments2 <- segregate(
  request = Input_tibble,
  RR = RecRates,
  MM = M_meta
)

Markers <- segments2markers(
  Segs = Segments2,
  Im = I_meta,
  Mm = M_meta,
  G = Geno
)

## -----------------------------------------------------------------------------
dim(Markers$ret_geno)

## -----------------------------------------------------------------------------
# genotypes of the first 10 individuals at the first 3 markers
Markers$ret_geno[1:10, 1:6]

## -----------------------------------------------------------------------------
# Individual IDs for the first 10 individuals
Markers$ret_ids[1:10,]

## -----------------------------------------------------------------------------
Input_tibble <- tibble(
  gpp = list(GSP, gsp4),
  reppop = list(RepPop1, RepPop4)
)

# here is what that input object looks like:
Input_tibble

## -----------------------------------------------------------------------------
set.seed(15) # for reproducibility
Segments <- segregate(
  request = Input_tibble,
  RR = RecRates
)

## ---- fig.width=7.5, fig.height=9---------------------------------------------
g <- plot_simulated_chromomsome_segments(Segments, RecRates)
g

## -----------------------------------------------------------------------------
Markers <- segments2markers(
  Segs = Segments,
  Im = I_meta,
  Mm = M_meta,
  G = Geno
)

## -----------------------------------------------------------------------------
####### These lines are not portable, since gscrambleTutorial.ped and
####### gscrambleTutorial.map are not available 
# with .ped / .map files

#plinkIN <- plink2gscramble(ped = "gscrambleTutorial.ped", map = "gscrambleTutorial.map")
#ls(plinkIN)
#str(plinkIN)

#dim(plinkIN$Geno)
#plinkIN$Geno[1:3,1:10]

#rm(plinkIN)

# identifical result with corresponding .ped.gz / .map.gz files
###### End non-portable lines

map_plink <- system.file("extdata/example-plink.map.gz", package = "gscramble")
ped_plink <- system.file("extdata/example-plink.ped.gz", package = "gscramble")

plinkIN <- plink2gscramble(ped_plink, map_plink)
dim(plinkIN$Geno)
dim(plinkIN$I_meta)
dim(plinkIN$M_meta)


## -----------------------------------------------------------------------------
prefix <- system.file("extdata/example-plink.map.gz", package = "gscramble")
prefix <- str_replace(prefix, "\\.map\\.gz", "")

plinkPRE <- plink2gscramble(prefix = prefix, gz_ext = TRUE)
dim(plinkPRE$Geno)
dim(plinkPRE$I_meta)
dim(plinkPRE$M_meta)

## -----------------------------------------------------------------------------
# note, we are only using a temporary file here for the output
# because this is in the vignette of an R package.  You, yourself,
# will likely want to write the results to a directory in your
# home directory somewhere, like `prefix = ~/my_stuff/gscram-sim-1`, etc.
tfile <- tempfile()
gscramble2plink(
  I_meta = Markers$ret_ids,
  M_meta = M_meta,
  Geno = Markers$ret_geno,
  prefix = tfile
)

