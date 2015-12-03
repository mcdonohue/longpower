# prepare the package for release

NEWS     = NEWS
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: news check clean

deps:
	Rscript -e 'if (!require("Rd2roxygen")) install.packages("Rd2roxygen", repos="http://cran.r-project.org")'

docs:
	R -q -e 'library(Rd2roxygen); rab(".", build = FALSE)'

roxy:
	R -q -e 'library(roxygen2); roxygenize(".")'

build: roxy
	cd ..;\
	R CMD build $(PKGNAME)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/;\
        $(RM) $(PKGNAME)_$(PKGVERS).tar.gz