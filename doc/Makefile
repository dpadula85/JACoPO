SHELL := /bin/bash

.PHONY: all view
all:
	pdflatex JACoPO.tex &> /dev/null
	bibtex JACoPO.aux &> /dev/null
	pdflatex JACoPO.tex &> /dev/null
	pdflatex JACoPO.tex &> /dev/null

view:
	xdg-open JACoPO.pdf

clean: 
	rm -rf *.out *.log *.aux *.bbl *.backup *.blg rsc*
