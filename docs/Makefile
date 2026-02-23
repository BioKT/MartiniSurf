SPHINXBUILD ?= sphinx-build
SOURCEDIR   = .
BUILDDIR    = _build

.PHONY: html clean

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html
	@echo "Built docs at $(BUILDDIR)/html/index.html"

clean:
	rm -rf $(BUILDDIR)
