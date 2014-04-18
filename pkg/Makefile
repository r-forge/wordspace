PKG_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' wordspace/DESCRIPTION)
PKG := wordspace_$(PKG_VERSION).tar.gz

INSTALL_FLAGS = 

# you may create a local file "local_config.mk" to override INSTALL_FLAGS or other settings
-include local_config.mk

default:
	@echo "Available targets for wordspace v$(PKG_VERSION)"
	@echo "  make build    ...  build source package"
	@echo "  make check    ...  test source package"
	@echo "  make install  ...  install source package locally"
	@echo "  make release  ...  check out and test source package for CRAN release"
	@echo "  make binary   ...  compile binary package for local architecture"
	@echo "  make clean    ...  clean up generated files"

.PHONY:	build check install release clean

build:
	Rscript --default-packages=Rcpp -e 'compileAttributes("wordspace")'
	R CMD build wordspace

check:	clean build $(PKG)
	R CMD check --multiarch $(PKG)

install:	clean build $(PKG)
	R CMD INSTALL $(INSTALL_FLAGS) $(PKG)

release:	clean
	R CMD build --resave-data wordspace
	R CMD check --as-cran $(PKG)
	@echo "$(PKG) is ready for release to CRAN"

binary:	build $(PKG)
	R CMD INSTALL --build $(PKG)

clean:
	-rm -f wordspace_*.tar.gz  # source package(s)
	-rm -f wordspace_*.tgz     # Mac binaries
	-rm -f wordspace_*.zip     # Windows binaries
	-rm -rf wordspace.Rcheck   # working directory of R CMD check
