## Updates
* June 14, 2015 - Version 1.3.0 - Major bug fix: version 1.2.0 did not fully resolve the issue with tapply(). Replaced with custom function.
* June 12, 2015 - Version 1.2.0 - Major bug fix: changed tapply() to ave() to preserve ordering of factors

## Test environments
* Local Windows 8 install, R 3.2.0
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs. 
* There were 2 NOTES on the local Windows check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Days since last update: 1
            * EXPLANATION: The previous minor version update did not fully resolve the bug. I believe the problem is now resolved. My apologies for not catching this in the most recent update!
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
* There were 2 NOTES on the win-builder (devel) check
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Days since last update: 1
            * EXPLANATION: The previous minor version update did not fully resolve the bug. I believe the problem is now resolved. 
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        phylogenetic (8:169)
            * EXPLANATION: This spelling is correct
* There was 1 NOTE on the win-builder (release) check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        * Days since last update: 1
            * EXPLANATION: The previous minor version update did not fully resolve the bug. I believe the problem is now resolved. 
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (4:8)
        phylogenetic (8:169)
            * EXPLANATION: This spelling is correct