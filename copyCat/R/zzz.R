.onLoad <- function(libname, pkgname) {
  if( !require(methods) ) stop("we require methods for package copyCat")
  initRdClass()
  print("Using copyCat version 1.5.3")
}
