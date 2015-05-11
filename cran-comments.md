## Test environments
* Local Windows 8 install, R 3.2.0
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs. 
* There was 1 NOTE on the local Windows check
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
* There were 2 NOTES on the win-builder (devel) check
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        phylogenetic (8:169)
            * EXPLANATION: This spelling is correct
* There was 1 NOTE on the win-builder (release) check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        phylogenetic (8:169)
            * EXPLANATION: This spelling is correct