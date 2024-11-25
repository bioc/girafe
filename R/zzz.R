
.onLoad <- function(libname, pkgname) {
  ## nothing to do here for the moment
}

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.22")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
  ## show vignette in windows menu
  if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("girafe")
  }
} # .onAttach

.onUnload <- function( libpath ) {
  ## nothing to do here for the moment
} # .onUnload
