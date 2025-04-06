# .onLoad <- function(libname, pkgname) {
#   # Re-assign the namespace of the defined functions to mgcv
#   assignInNamespace("fix.family.link.family", fix.family.link.family, ns = "mgcv")
#   assignInNamespace("fix.family.var", fix.family.var, ns = "mgcv")
#
# }
.onLoad <- function(libname, pkgname) {
  # Unlock and replace locked bindings
  unlockBinding("fix.family.link.family", asNamespace("mgcv"))
  assign("fix.family.link.family", fix.family.link.family, envir = asNamespace("mgcv"))
  lockBinding("fix.family.link.family", asNamespace("mgcv"))

  unlockBinding("fix.family.var", asNamespace("mgcv"))
  assign("fix.family.var", fix.family.var, envir = asNamespace("mgcv"))
  lockBinding("fix.family.var", asNamespace("mgcv"))
}
