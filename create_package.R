# remove old version first
remove.packages('nmmsims')

# create R functions from C++ files
Rcpp::compileAttributes()
# nmmsims_package_exists <- ("nmmsims" %in% installed.packages()[,"Package"])

# Say that the package uses UTF-8 for stuff
# write('Encoding: UTF-8', file='DESCRIPTION', append=TRUE)

# create documentation
devtools::document()

devtools::build(path = 'G:/My Drive/Sharing')

# install the package
devtools::install(dependencies=TRUE)

library('nmmsims')
# nmmsims::

# download.file(url = 'https://cdn-33.anonfile.com/zaOeO4Dfn8/8cbe1a96-1575689524/nmmsims_1.0.tar.gz', destfile = 'nmmsims_1.0.tar.gz')
# install.packages('nmmsims_1.0.tar.gz', repos = NULL, type="source")
# install.packages('https://cdn-33.anonfile.com/zaOeO4Dfn8/8cbe1a96-1575689524/nmmsims_1.0.tar.gz', repos = NULL, type="source")

# there are many more potential improvements here, see: http://r-pkgs.had.co.nz/src.html
# Whenever you use C++ code in your package, you need to clean up after yourself when your package is unloaded. 
# Do this by writing a .onUnload() function that unloads the DLL:
.onUnload <- function (libpath) {
  library.dynam.unload("nmmsims", libpath)
}
# but... https://stackoverflow.com/questions/26691878/must-r-packages-unload-dynamic-libraries-when-they-unload

# devtools::install_github(repo = 'philipshirk/nmmsims')
