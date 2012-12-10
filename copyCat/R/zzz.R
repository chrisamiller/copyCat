.onLoad <- function(libname, pkgname) {
  if( !require(methods) ) stop("we require methods for package readDepth")
  initRdClass()
  registerDoMC()
  print("Using readDepth version 1.5.0")
}
