.onLoad <- function(libname, pkgname) {
  # register custom gating functions on load so ready to be used internal to openCyto
  openCyto::registerPlugins(fun=.specificMindensity,methodName='specificMindensity',dep=NA)
  openCyto::registerPlugins(fun=.dummySubgate,methodName='dummySubgate',dep=NA)
  
}