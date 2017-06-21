# TODO: Add comment
# 
# Author: klambaue
###############################################################################


.onLoad <- function(libname, pkgname) {
	library.dynam("gcn.mops", pkgname, libname)
}

.onUnload <- function(libpath)
{
	library.dynam.unload("gcn.mops", libpath)
}
