.onLoad <- function(libname, pkgname) {
  if( !require(methods) ) stop("we require methods for package readDepth")
  initRdClass()
  registerDoMC()
  print("Using readDepth version 1.3.1")
}
