SHELL=/bin/bash
NOMBRE=$(shell grep -H \\documentclass *.tex | cut -d: -f1 | cut -d. -f1)

ifeq ($(strip $(wildcard $(NOMBRE).bib) ),)
HAVE_BIB=no
else
HAVE_BIB=yes
endif

DEPENDENCIES=$(NOMBRE).tex

ifeq '$(HAVE_BIB)' 'yes'
DEPENDENCIES+= $(NOMBRE).bib
endif

.PHONY: clean distclean all k 

all: $(NOMBRE).pdf

$(NOMBRE).pdf: $(DEPENDENCIES)
	@echo "   Making dvi for first time..."
	@latex -halt-on-error $(NOMBRE).tex &> .my_log || (cat .my_log && rm .my_log && exit 1)
	@rm .my_log
	@if [[ "${HAVE_BIB}" == "yes" ]]; then \
	 echo "   Compiling references..." && \
	 (bibtex $(NOMBRE) &> .my_log || (cat .my_log && rm .my_log && exit 1)) && \
	 rm .my_log && \
	 echo "   Re-making dvi including bibliography..." && \
	 (latex -halt-on-error $(NOMBRE).tex &> .my_log && rm .my_log) || (cat .my_log && rm .my_log && exit 1) \
	fi
	@echo "   Re-making dvi for satisfying references..."
	@latex $(NOMBRE).tex &> /dev/null
	@echo "   Generating final pdf..."
	@dvipdf $(NOMBRE).dvi &> /dev/null

clean:
	-rm -f $(NOMBRE).{aux,toc,log,tmp,dvi,idx,ilg,ind,bbl,blg,out} .my_log

distclean: clean
	-rm -f $(NOMBRE).pdf

k: $(NOMBRE).pdf
	@echo "   Opening $(NOMBRE).pdf with kpdf..."
	@kpdf $(NOMBRE).pdf &> /dev/null &

x: $(NOMBRE).pdf
	@echo "   Opening $(NOMBRE).pdf with xpdf..."
	@evince $(NOMBRE).pdf &> /dev/null &

e: $(NOMBRE).pdf
	@echo "   Opening $(NOMBRE).pdf with evince..."
	@evince $(NOMBRE).pdf &> /dev/null &
