## Test environments
* Local Windows 10 install, R 3.3.2
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs, WARNINGs, or NOTEs.
* There was 1 NOTE on the win-builder (devel and release) check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        morphometric (7:75)
        phylogenetic (7:165)
            * EXPLANATION: This spelling is correct
        * Uses the superseded package: 'doSNOW'
            * EXPLANATION: imports doSNOW and doParallel for cross-platform compatability