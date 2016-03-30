PKG_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' wordspace/DESCRIPTION)
PKG := wordspace_$(PKG_VERSION).tar.gz

SSVD_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' sparsesvd/DESCRIPTION)
SSVD := sparsesvd_$(SSVD_VERSION).tar.gz

PKGS = $(PKG) $(SSVD)

INSTALL_FLAGS = 

# you may create a local file "local_config.mk" to override INSTALL_FLAGS or other settings
-include local_config.mk

default:
	@echo "Available targets:"
	@echo "  make build[_pkg]    ...  build source package(s)"
	@echo "  make check[_pkg]    ...  test source package(s)"
	@echo "  make install[_pkg]  ...  install source package locally(s)"
	@echo "  make release[_pkg]  ...  check out and test source package(s) for CRAN release"
	@echo "  make binary[_pkg]   ...  compile binary package(s) for local architecture"
	@echo "  make clean[_pkg]    ...  clean up generated files"
	@echo "Packages: wordspace (v$(PKG_VERSION)), sparsesvd (v$(SSVD_VERSION))"

.PHONY:	build check install release clean
.PHONY:	build_wordspace check_wordspace install_wordspace release_wordspace clean_wordspace
.PHONY:	build_sparsesvd check_sparsesvd install_sparsesvd release_sparsesvd clean_sparsesvd


build:	build_wordspace build_sparsesvd

build_wordspace:
	cat wordspace/src/parts/HEADER wordspace/src/parts/*.h wordspace/src/parts/*.cpp > wordspace/src/wordspace.cpp
	Rscript --default-packages=Rcpp -e 'compileAttributes("wordspace")'
	R CMD build wordspace

build_sparsesvd:
	R CMD build sparsesvd


check:	check_wordspace check_sparsesvd

check_wordspace:	clean_wordspace build_wordspace $(PKG)
	R CMD check --multiarch $(PKG)

check_sparsesvd:	clean_sparsesvd build_sparsesvd $(SSVD)
	R CMD check --multiarch $(SSVD)


install:	install_wordspace install_sparsesvd

install_wordspace:	clean_wordspace build_wordspace $(PKG)
	R CMD INSTALL $(INSTALL_FLAGS) $(PKG)

install_sparsesvd:	clean_sparsesvd build_sparsesvd $(SSVD)
	R CMD INSTALL $(INSTALL_FLAGS) $(SSVD)


release:	release_wordspace release_sparsesvd

release_wordspace:	clean_wordspace
	R CMD build --resave-data wordspace
	R CMD check --as-cran $(PKG)
	@echo "$(PKG) is ready for release to CRAN"

release_sparsesvd:	clean_sparsesvd
	R CMD build --resave-data sparsesvd
	R CMD check --as-cran $(SSVD)
	@echo "$(SSVD) is ready for release to CRAN"


binary:	binary_wordspace binary_sparsesvd

binary_wordspace:	build_wordspace $(PKG)
	R CMD INSTALL --build $(PKG)

binary_sparsesvd:	build_sparsesvd $(SSVD)
	R CMD INSTALL --build $(SSVD)


clean:	clean_wordspace clean_sparsesvd

clean_wordspace:
	-rm -f  wordspace_*.tar.gz  # source package(s)
	-rm -f  wordspace_*.tgz     # Mac binaries
	-rm -f  wordspace_*.zip     # Windows binaries
	-rm -rf wordspace.Rcheck    # working directory of R CMD check

clean_sparsesvd:
	-rm -f  sparsesvd_*.tar.gz  # source package(s)
	-rm -f  sparsesvd_*.tgz     # Mac binaries
	-rm -f  sparsesvd_*.zip     # Windows binaries
	-rm -rf sparsesvd.Rcheck    # working directory of R CMD check
