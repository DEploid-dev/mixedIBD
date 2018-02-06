.PHONY: all clean
all: clean main.pdf

mainfiguresPng = $(shell grep png main.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
mainfiguresPdf = $(shell grep pdf main.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
supfigures = $(shell grep png mainSupplement.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
supfigurespdf = $(shell grep pdf mainSupplement.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
suptex = $(shell grep "\.tex" mainSupplement.tex | sed -e "s/^.*{//g" -e "s/\}//g" )

coverLetter.pdf: coverLetter.tex
	pdflatex coverLetter.tex

main.pdf: main.tex ${mainfiguresPng} ${mainfiguresPdf} main.bbl
	pdflatex main.tex
	pdflatex main.tex

main.aux: main.tex
	pdflatex main.tex

main.bbl: main.aux
	bibtex main.aux

mainSupplement.pdf: mainSupplement.tex ${supfigures} ${supfigurespdf} ${suptex} supplementReset.tex
	pdflatex mainSupplement.tex
	pdflatex mainSupplement.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out main.pdf

