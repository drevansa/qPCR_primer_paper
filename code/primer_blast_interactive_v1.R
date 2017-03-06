#!/usr/bin/env Rscript

################################################################################
## primer_blast_v1.R
################################################################################
## Description: Batch processing primer pair submission to NCBI Primer-BLAST
##==============================================================================
## Author: A. Evans 22/02/2017
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## locale: en_GB.UTF-8
##------------------------------------------------------------------------------
##==============================================================================
## Requires input CSV file with format:
##==============================================================================
## NAME,PAIR,FWRD,REV
## Source 1,Pair1,GTGCCAGCAGCCGCGGTAA,GGACTACAAGGGTATCTAAT
## Source 1,Pair2,GTGCCAGCCGCCGCGGTAA,GGACTACAAGGGTATCTAAT
## Source 2,Pair1,GTGTCAGCAGCCGCGGTAA,GGACTACAAGGGTATCTAAT
## 
## Where: 
## Header requires columns named: NAME, PAIR, FWRD and REV
## Name     - Text string
## PAIR     - Text string
## FWRD     - Forward primer
## REV      - Reverse primer
##
## Additional columns are permitted, but will not be imported.
##
##==============================================================================


#################################################
## SET VARIABLES
#################################################
#setwd("PATH")
#infl <- "PATH/primer_blast_inputs.csv"
#outfl <- "PATH/primer_blast_outputs.csv"
#################################################


## library check and load

deps <- c("methods", "optparse", "primerTree", "plyr", "seqinr", "directlabels",
          "gridExtra")

missingdeps <- deps[!(deps %in% installed.packages()[ ,"Package"])]

if(length(missingdeps)) {
    missingdeps <- paste(missingdeps,collapse="; ")
    stop(paste("Install packages:", missingdeps, sep = " "))
}


library("methods", quietly = TRUE)
library("optparse", quietly = TRUE)
library("primerTree", quietly = TRUE)
library("plyr", quietly = TRUE)
library("seqinr", warn.conflicts = FALSE, quietly = TRUE)


## Increasing limit - else search results are truncated

options(max.print = 1000000)


## If input file exists, build data frame of primer pairs

if (is.null(infl)){
    print_help(opt_parser)
    print("Supply input csv file name [Rscript primer_blast_v1.R -i filename]")
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    stop()

} else if (!file.exists(infl)) {
    msg <- paste("File does not exist: ", infl, sep="")
    stop(msg)

} else {
    print("Script called with the following:")
    print(paste("input: ", infl, sep = ""))
    print(paste("output: ", outfl, sep = ""))

    nr <- length(count.fields(infl, sep = "\n"))

    vars <- rep("NULL", 4)
    vars[1] <- "character"   # NAME
    vars[2] <- "character"   # PAIR
    vars[3] <- "character"   # FWRD
    vars[4] <- "character"   # REV

    cc = vars

    primerpairs <- read.csv(infl, sep = ",", na.strings = "", header = TRUE,
                            nrows = nr, colClasses = cc, strip.white = TRUE)
}


## Set header for output file

hdr<-"SheetSort,Hit,Source,Pair,FWD_Primer,REV_Primer,FWD_Primer_Len,REV_Primer_Len,FWD_Primer_GC,REV_Primer_GC,gi,accession,species,product_length,mismatch_forward,mismatch_reverse,forward_start,forward_stop,reverse_start,reverse_stop,product_start,product_stop,Date,NCBI_URL"

writeLines(hdr, outfl)


## Sumbit primer pairs to Primer-Blast
## NOTE: SALT_FORMULAR="1" and TM_METHOD="1" >  SantaLucia 1998
##       We are restricting num_aligns to 25
##       Submitting one job at a time to avoid server overload
##       This may take a while!

si <- 1

totrows <- 0
totpp <- nrow(primerpairs)


for (i in 1:(nr -1)) {

    IN_NAME <- paste(primerpairs$NAME[i], primerpairs$PAIR[i], sep = ",")
    IN_PRIMER_F <- primerpairs$FWRD[i]
    IN_PRIMER_R <- primerpairs$REV[i]

    print(paste("Processing:", i, "of", totpp, IN_NAME, IN_PRIMER_F, IN_PRIMER_R, sep=" "))

    SEARCH_RES = search_primer_pair(name = IN_NAME,
    simplify = TRUE,
    IN_PRIMER_F,
    IN_PRIMER_R,
    num_aligns = 25,
    PRIMER_PRODUCT_MIN = "70",
    PRIMER_PRODUCT_MAX = "1000",
    PRIMER_NUM_RETURN = "10",
    PRIMER_MIN_TM = "57.0",
    PRIMER_OPT_TM = "60.0",
    PRIMER_MAX_TM = "63.0",
    PRIMER_MAX_DIFF_TM = "3",
    SPLICE_SITE_OVERLAP_5END = "7",
    SPLICE_SITE_OVERLAP_3END = "4",
    MIN_INTRON_SIZE = "1000",
    MAX_INTRON_SIZE = "1000000",
    primer_specificity_database = "refseq_representative_genomes",
    exclude_xm = "off",
    exclude_env = "off",
    organism = "Fungi",
    TOTAL_PRIMER_SPECIFICITY_MISMATCH = "2",
    PRIMER_3END_SPECIFICITY_MISMATCH = "2",
    MISMATCH_REGION_LENGTH = "5",
    TOTAL_MISMATCH_IGNORE = "6",
    MAX_TARGET_SIZE = "4000",
    HITSIZE = "50000",
    EVALUE = "100000",
    WORD_SIZE = "6",
    MAX_CANDIDATE_PRIMER = "500",
    NUM_TARGETS = "20",
    MAX_TARGET_PER_TEMPLATE = "1000",
    PRIMER_MIN_SIZE = "15",
    PRIMER_OPT_SIZE = "20",
    PRIMER_MAX_SIZE = "25",
    PRIMER_MIN_GC = "20.0",
    PRIMER_MAX_GC = "80.0",
    GC_CLAMP = "0",
    POLYX = "5",
    PRIMER_MAX_END_STABILITY = "9",
    PRIMER_MAX_END_GC = "5",
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH = "40.00",
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH = "70.00",
    PRIMER_MAX_SELF_ANY_TH = "45.0",
    PRIMER_MAX_SELF_END_TH = "35.0",
    PRIMER_PAIR_MAX_COMPL_ANY_TH = "45.0",
    PRIMER_PAIR_MAX_COMPL_END_TH = "35.0",
    PRIMER_MAX_HAIRPIN_TH = "24.0",
    PRIMER_MAX_TEMPLATE_MISPRIMING = "12.00",
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING = "24.00",
    SELF_ANY = "8.00",
    SELF_END = "3.00",
    PRIMER_PAIR_MAX_COMPL_ANY = "8.00",
    PRIMER_PAIR_MAX_COMPL_END = "3.00",
    OVERLAP_5END = "7",
    OVERLAP_3END = "4",
    MONO_CATIONS = "50.0",
    DIVA_CATIONS = "1.5",
    CON_DNTPS = "0.6",
    SALT_FORMULAR = "1",
    TM_METHOD = "1",
    CON_ANEAL_OLIGO = "50.0",
    PRIMER_MISPRIMING_LIBRARY = "AUTO",
    PRIMER_INTERNAL_OLIGO_MIN_SIZE = "18",
    PRIMER_INTERNAL_OLIGO_OPT_SIZE = "20",
    PRIMER_INTERNAL_OLIGO_MAX_SIZE = "27",
    PRIMER_INTERNAL_OLIGO_MIN_TM = "57.0",
    PRIMER_INTERNAL_OLIGO_OPT_TM = "60.0",
    PRIMER_INTERNAL_OLIGO_MAX_TM = "63.0",
    PRIMER_INTERNAL_OLIGO_MIN_GC = "20.0",
    PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT = "50",
    PRIMER_INTERNAL_OLIGO_MAX_GC = "80.0" )


    ## If Blast Fails exit, else append results to outfile

    if (is.null(SEARCH_RES$BLAST_result)){
        msg <- (paste("Failed to retrieve Primer-Blast results: ", i, "of", totpp, 
                       IN_NAME, IN_PRIMER_F, IN_PRIMER_R, sep = " "))
        opt <- options(show.error.messages = TRUE)
        on.exit(options(opt))
        stop(msg)
    }  else if (length(SEARCH_RES$taxonomy$gi) < 1) {
        msg <- (paste("Failed to retrieve taxonomy GI: ", i, "of", totpp, IN_NAME, 
                       IN_PRIMER_F, IN_PRIMER_R, sep = " "))
        opt <- options(show.error.messages = TRUE)
        on.exit(options(opt))
        stop(msg)
    } else if (is.null(SEARCH_RES$taxonomy$species)) {
        msg <- (paste("Failed to retrieve species: ", i, "of", totpp, IN_NAME, 
                       IN_PRIMER_F, IN_PRIMER_R, sep = " "))
        opt <- options(show.error.messages = TRUE)
        on.exit(options(opt))
        stop(msg)
    } else {

        ## Extract results and append to outfile

        df1<-data.frame(Source = primerpairs$NAME[i],
                        Pair = primerpairs$PAIR[i],
                        FWD_Primer = IN_PRIMER_F,
                        REV_Primer = IN_PRIMER_R,
                        FWD_Primer_Len = nchar(IN_PRIMER_F),
                        REV_Primer_Len = nchar(IN_PRIMER_R),
                        FWD_Primer_GC = GC(unlist(strsplit(IN_PRIMER_F,split = NULL))),
                        REV_Primer_GC = GC(unlist(strsplit(IN_PRIMER_R,split = NULL))),
                        NCBI_URL = SEARCH_RES$response$`1`$url,
                        Date = SEARCH_RES$response$`1`$date)

        df2 <- data.frame(SEARCH_RES$BLAST_result)

        df3 <- data.frame(gi = SEARCH_RES$taxonomy$gi,species = SEARCH_RES$taxonomy$species)

        df4 <- join(df2, df3, by = 'gi', type = 'left', match = 'all')

        totrows <- (si + nrow(df4)) -1
        ssdf <- data.frame(sheetsort = seq(si, totrows, by = 1))
        hitdf<-data.frame(hit = seq(1, nrow(df4), by = 1))

        df5 <- cbind(ssdf, hitdf, df1, df4)

        dfout <- subset(df5, select = c(sheetsort,
                                      hit,
                                      Source,
                                      Pair,
                                      FWD_Primer,
                                      REV_Primer,
                                      FWD_Primer_Len,
                                      REV_Primer_Len,
                                      FWD_Primer_GC,
                                      REV_Primer_GC,
                                      gi,
                                      accession,
                                      species,
                                      product_length:product_stop,
                                      Date,
                                      NCBI_URL)
                       )

        write.table(dfout, file = outfl, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

        # increment sheet sort index
        si <- si + nrow(df4)
    }
    rm(SEARCH_RES)
}

