all:	pdf

pdf:	girafe.tex
	pdflatex girafe
	bibtex girafe
	pdflatex girafe
	pdflatex girafe
	make clean

clean:
	rm -f *.out *.bbl *.aux *.log *.blg *.brf *.toc  *~ *.eps
	rm -f *.nfs*
	rm -f girafe-*.png girafe-*.pdf
#	rm -f Rplots.ps Rplots.pdf   # problem on build machine gewurz
#	rm -f *.dvi *.ps girafe.tex
