PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
BIOCVER := RELEASE_3_23
QUARTO_CACHE := $(shell pwd)/.quarto-cache
QUARTO_LIB := $(shell pwd)/.quarto-lib
QUARTO_R_LIBS := $(shell Rscript -e 'cat(paste(c(normalizePath(".quarto-lib", winslash = "/", mustWork = FALSE), .libPaths()), collapse = .Platform$$path.sep))')
PUBLISH_DIR ?= .site


all: rd check clean

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

vignette:
	mkdir -p $(QUARTO_CACHE);\
	mkdir -p $(QUARTO_LIB);\
	mkdir -p $(PUBLISH_DIR);\
	R CMD INSTALL -l $(QUARTO_LIB) .;\
	cd vignettes;\
	R_LIBS="$(QUARTO_R_LIBS)" XDG_CACHE_HOME=$(QUARTO_CACHE) LOCALAPPDATA=$(QUARTO_CACHE) quarto render seqcombo.qmd --to html;\
	mv seqcombo.html ../$(PUBLISH_DIR)/index.html

build:
	cd ..;\
	R CMD build $(PKGSRC)

build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

rmrelease:
	git branch -D $(BIOCVER)

release:
	git checkout $(BIOCVER);\
	git fetch --all

update:
	git fetch --all;\
	git checkout devel;\
	git merge upstream/devel;\
	git merge origin/devel

push: update
	git push upstream devel;\
	git push origin devel

biocinit:
	git remote add upstream git@git.bioconductor.org:packages/$(PKGNAME).git;\
	git fetch --all


