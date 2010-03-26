
.onLoad <- function(libname, pkgname) {
  ## nothing to do here for the moment
}

.onAttach <- function(libname, pkgname) {
  ## show vignette in windows menu
  if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("girafe")
  }
} # .onAttach

.onUnload <- function( libpath ) {
  ## nothing to do here for the moment
} # .onUnload
