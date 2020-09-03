## Resubmission
This is a resubmission. In this version I have:
* Fixed redundancies in DESCRIPTION
* Added a doi reference to the literature

## Test environments
* local MacOS Cadtalina install, (R 4.0.2)
* ubuntu 16.04 (Github R actions), (R release, old-release and 3.5)
* win-builder (devel, release and old-release)
* Using R-hub:-  

      *  Debian Linux, R-release, GCC  
      *  Debian Linux, R-devel, clang, ISO-8859-15 locale  
      *  Debian Linux, R-devel, GCC  
      *  Fedora Linux, R-devel, clang, gfortran  
      *  Fedora Linux, R-devel, GCC  
      *  macOS 10.13.6 High Sierra, R-release, CRAN's setup  

## R CMD check results
There were no ERRORs or WARNINGs. 

Here, I have ignored the NOTE on "new submission".

There was 1 NOTE:

* checking for future file timestamps ... NOTE  
  unable to verify current time.

I understand that this is an issue with the unavailability of <http://worldclockapi.com/> that R CMD check uses, and is not relevant.