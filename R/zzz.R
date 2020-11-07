.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Dune now uses the Normalized Mutual Information",
    "instead of the adjusted Rand Index. You can restore the old option by",
    "setting metric'ARI' in the Dune function."))
}