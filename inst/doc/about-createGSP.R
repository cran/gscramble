## ----setup, include = FALSE---------------------------------------------------
library(gscramble)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# determine if dot is on the system or not.  If we are running it on
# our own machines, with dot, then it will update the images. If dot is
# not on the system it will just use the stored dot images.
HasDot <- Sys.which("dot") != ""


# function to get the bits of the index
bits <- function(i) {as.logical(sapply(i,function(x){ as.logical(intToBits(x))}))[1:4]}

# a vector of the hybrid categories possible
classes <- c("F1", "F2", "F1B", "F1B2")

# a function to create (if dot is present) the figure for the i-th
# component of GSP_opts
make_plots <- function(i) {
  
  # get the section title for it
  section_title <- paste(
    "### ",
    #names(GSP_opts)[i],
    #":  ",
    paste(classes, " = ", bits(i), ", &nbsp;&nbsp;&nbsp;&nbsp; ", sep = "", collapse = "")
  )
  
  # make the plots if dot is installed
  if(HasDot) {
    gsp2dot(GSP_opts[[i]], paste("../man/figures/", i, sep=""))
  }
  
  # remove the dot and png files, that we do not need
  catch <- file.remove(paste("../man/figures/", i, c(".dot", ".png"), sep = ""))
  
  # now actually make the text for the sections
  ticks3 <- "```"
  sprintf("\n\n%s\n\n%s{r, echo=FALSE, out.width='100%%', fig.align='center'}\nknitr::include_graphics('../man/figures/%d.svg')\n%s", section_title, ticks3, i, ticks3)
}

# uncomment the following line and run in the console to get text for 
# all the sections:
#dump <- lapply(1:15, function(i) cat(make_plots(i)))

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/1.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/2.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/3.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/4.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/5.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/6.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/7.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/8.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/9.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/10.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/11.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/12.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/13.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/14.svg')

## ---- echo=FALSE, out.width='100%', fig.align='center'------------------------
knitr::include_graphics('../man/figures/15.svg')

