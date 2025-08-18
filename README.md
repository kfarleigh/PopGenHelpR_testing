# PopGenHelpR_testing
This repository contains data and code to test the R package [PopGenHelpR](https://cran.r-project.org/web/packages/PopGenHelpR/index.html). This workflow could not be included with the package because it is so large. Details on files are below.

The testing will take hours if not days, to run yourself. You are welcome to, but we wanted to save you time. 

## Markdown report
The [PopGenHelpR_testing.html](https://github.com/kfarleigh/PopGenHelpR_testing/blob/main/PopGenHelpR_testing.html) is an html file that contains the output of the testing. This was rendered and knitted using the `PopGenHelpR_testing.Rmd` file. 

## Script

The `PGH_Stattesting_forMS.R` is the same as the .Rmd file, but it also contains tests to write output files and ensure the `Private.alleles` function works with phased and unphased data. 

## Rdata file
The `PGH_stattesting_workspace.RData` contains the output of `PGH_Stattesting_forMS.R`  in case you would like to check the results for yourself. We supply this because the testing can take hours to days depending on your machine. 

Please reach out to Keaka Farleigh if you have any questions (see email on CRAN page linked above). 
