## Resubmission

The package was accepted to CRAN on 10th Sep 2020, URL :
  https://cran.r-project.org/web/packages/qbld/index.html.


However, an email was received almost immediately that there was 
an error in the CRAN check (only on Solaris OS), and I needed to resubmit by 24th.

In this version update (patch), I have identified and fixed issues as per 
the install log and the package build has been checked on r-hub additionally 
on the following platforms :-

- Oracle Solaris 10, x86, 32 bit,
  R-release
- Oracle Solaris 10, x86, 32 bit,
  R-release, Oracle Developer Studio 12.6

No ERROR(s), WARNING(s) OR NOTE(s) found.



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