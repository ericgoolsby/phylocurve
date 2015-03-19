## Test environments
* Local Windows 8 install, R 3.1.2
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs. 
* There were 2 NOTES on the local Windows check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
          New submission
            * EXPLANATION: This is a new submission
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
* There was 1 NOTE on the win-builder (devel and release) check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        New submission
            * EXPLANATION: This is a new submission
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        phylogenetic (8:169)
            * EXPLANATION: This spelling is correct