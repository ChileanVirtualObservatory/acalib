default: O1-5.pdf 

clean:
	rm -f *~ *.aux *.bbl *.blg *.lof *.log *.lot *.toc *.nav *.snm *.idx *.out 

distclean: clean
	rm -f *.pdf *.dvi

.PHONY: all clean distclean

%.dvi: %.tex
	latex -halt-on-error $<
	bibtex $(basename $<)
	latex -halt-on-error $<
	latex -halt-on-error $<

%.pdf: %.dvi
	dvipdf $<
