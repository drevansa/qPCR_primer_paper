#!/usr/bin/env Rscript

################################################################################
## primer_blast_mismatch_plots_v1.R
################################################################################
## Description: Generate 3 x 3 primer pair mismatch frequency plots
##==============================================================================
## Author: A. Evans 24/02/2017
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## locale: en_GB.UTF-8
##------------------------------------------------------------------------------
##==============================================================================
## Input
##==============================================================================
## Takes primer_blast_v1.R output csv file [24 col by x rows]
## generated by study research criteria [9 x sources]
##
## Reads in:
##      col 3  Header: Source
##      col 13 Header: species
##      col 15 Header: mismatches_forward
##      col 16 Header: mismatches_reverse
##
##==============================================================================
## Output
##==============================================================================
## 3 x 3 bar plots: .pdf, .tiff and .eps
## mismatch frequency by source table:.txt
##==============================================================================


#################################################
## SET VARIABLES
#################################################
#setwd("PATH")
#infl <- "PATH/primer_blast_outputs.csv"
#outfl <- "PATH/> getwd()
#################################################

## append extensions for output

outfl_tiff <- paste(outfl, ".tiff", sep = "")
outfl_eps <- paste(outfl, ".eps", sep = "")
outfl_pdf <- paste(outfl, ".pdf", sep = "")
outfl_txt <- paste(outfl, ".txt", sep = "")



## Read in four cols only

hdr <- read.csv(file = infl, nrows = 1)
n_cols <- length(hdr)

vars <- rep("NULL", n_cols)
vars[3]  <- "character"   # source
vars[13] <- "character"   # species
vars[15] <- "integer"     # mismatches_forward
vars[16] <- "integer"     # mismatches_reverse

cc = vars

pbdf <- read.csv(file = infl, header = TRUE, sep = ",", na.strings = "", colClasses = cc)



## Summarize mismatch counts

pbdf$tot_mismatches <- rowSums(pbdf[ ,3:4])

mismatch_tbl <- with(pbdf, table(Source, tot_mismatches))



## Amend for plot headers

row.names(mismatch_tbl) <- gsub("Caporaso", "Caporaso et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Fierer", "Fierer et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Heuer", "Heuer et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("HMP", "HMP", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Maeda", "Maeda et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Muyzer", "Muyzer et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Schwieger", "Schwieger \\& Tebbe", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Shakya", "Shakya et al", row.names(mismatch_tbl))
row.names(mismatch_tbl) <- gsub("Weisburg", "Weisburg et al", row.names(mismatch_tbl))


## generate 3 x 3 plots

par(mfrow = c(3, 3), mar = c(2.8, 2.6, 2.5, 1.0))

for (i in 1:9) {
    mlabel <- row.names(mismatch_tbl)[i]
    ylimits <- range(pretty(c(0, mismatch_tbl[i, ])))
    xlabel <- "Number of mismatches"
    ylabel <- "Frequency"
    cols <- c("black", "white")[(mismatch_tbl[i, ] < 1) + 1]

    mismatch_plot<- barplot(mismatch_tbl[i, ],  main = mlabel, yaxs = "r", 
                            ylim = ylimits, xlab = xlabel, ylab = ylabel, las = 1, 
                            col = cols, border = cols, ps = 7, cex.main = 1, 
                            cex.lab = 1, cex.axis = 0.86, cex.names = 0.86, 
                            tck = -0.025,  mgp=c(1.5, .3, 0))

    axis(1, at = mismatch_plot, labels = FALSE, tick = TRUE, ps = 7,  tck = -0.025)

    box(lwd = 0.8)

}


## output plots

dev.print(pdf, outfl_pdf)
dev.print(cairo_ps, outfl_eps)
dev.print(tiff, outfl_tiff, res = 600, family = "Arial", width = 190, height = 160, units = "mm")

dev.off()


## output frequency table
write.ftable(ftable(mismatch_tbl), file = outfl_txt, sep = "\t", quote = FALSE)

