######################################################################
#
# zzz.R
#
# Written by Ayush Agarwal <ayush.agarwal50@gmail.com>.
#
# Part of the R/qbld package
#
# .onLoad is run when the package is loaded with library(qbld)
#
######################################################################

#' @import utils

.onAttach = function(libname, pkgname)
{
  temp = packageDescription("qbld")
  msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
              " created on ", temp$Date, ".\n", sep = "")
  msg = paste(msg, '\nFor citation information, type citation("qbld").\n',sep = "")
  msg = paste(msg, 'Type help("qbld-package") or help("model.qbld") to get started.\n',sep = "")
  packageStartupMessage(msg)
}
