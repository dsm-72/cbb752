#!/usr/bin/env Rscript

### Usage: Rscript --vanilla hw1.R <input file> <score file>
### Example: Rscript --vanilla hw1.R input.txt blosum62.txt
### Note: Smith-Waterman Algorithm

### This is one way to read in arguments in R
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
stop("At least two arguments must be supplied (inputFile, scoreFile).n", call.=FALSE) } else if (length(args)>=2) {
    # default gap penalties
    args[3] = -2
    args[4] = -1 }

## Specifying author and email
p <- c(person("FirstName", "LastName", role = "aut", email = "first.last@yale.edu"))

## Implement your Smith-Waterman Algorithm
runSW <- function(inputFile, scoreFile, openGap = -2, extGap = -1) {
    ### calculation
    ### write output
}

## Run the main function and generate results
runSW(inputFile=args[1], scoreFile=args[2], openGap=args[3], extGap=args[4])