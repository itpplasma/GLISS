# Build the HTML documentation: make -C docs html
SPHINXBUILD ?= sphinx-build

html:
	$(SPHINXBUILD) -W -b html . _build/html

linkcheck:
	$(SPHINXBUILD) -b linkcheck . _build/linkcheck

clean:
	rm -rf _build

.PHONY: html linkcheck clean
