.onLoad <- function(libname, pkgname) 
{
  ##library.dynam(file.path(libname, paste("eaf", .Platform$dynlib.ext, sep="")), pkgname, libname )
  library.dynam("eaf", pkgname, libname )
}

#.noGenerics <- TRUE

.onUnload <- function(libpath)
  library.dynam.unload("eaf", libpath)




## If present, .First.lib will be used if the NAMESPACE file is
## missing.  This is useful during development, thanks to C-c C-l in
## Emacs/ESS. It won't be used if a NAMESPACE file is present. 
.First.lib <- function(lib, pkg) 
{
  library.dynam("eaf", pkg, lib )
}





