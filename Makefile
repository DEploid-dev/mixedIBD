.PHONY: all clean pdfFigures
all: clean main.pdf

mainfiguresPng = $(shell grep png main.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
mainfiguresPdf = $(shell grep "\.pdf" main.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
#supfigures = $(shell grep png mainSupplement.tex | sed -e "s/^.*{/figures\/otherFigures\//g" -e "s/\}//g" )
#supfigurespdf = $(shell grep pdf mainSupplement.tex | sed -e "s/^.*{/figures\/otherFigures\//g" -e "s/\}//g" )
suptex = $(shell grep "\.tex" mainSupplement.tex | sed -e "s/^.*{//g" -e "s/\}//g" )

coverLetter.pdf: coverLetter.tex
	pdflatex coverLetter.tex

main.pdf: main.tex ${mainfiguresPng} ${mainfiguresPdf} main.bbl pdfFigures
	pdflatex main.tex
	pdflatex main.tex

main_todo.pdf: main.pdf
	sed -e "s/\\\usepackage\[disable\]{todonotes}/\\\usepackage\[colorinlistoftodos\]{todonotes}/" \
	 -e "s/\\\textcolor{black}/\\\textcolor{red}/g" main.tex > main_todo.tex
	pdflatex main_todo.tex
	bibtex main_todo.aux
	pdflatex main_todo.tex
	pdflatex main_todo.tex

main.aux: main.tex
	pdflatex main.tex

main.bbl: main.aux
	bibtex main.aux

mainSupplement.pdf: mainSupplement.tex ${supfigures} ${supfigurespdf} ${suptex} supplementReset.tex mainSupplement.bbl
	pdflatex mainSupplement.tex
	pdflatex mainSupplement.tex

mainSupplement_todo.pdf: mainSupplement.pdf
	sed -e "s/\\\usepackage\[disable\]{todonotes}/\\\usepackage\[colorinlistoftodos\]{todonotes}/" \
	 -e "s/\\\textcolor{black}/\\\textcolor{red}/" mainSupplement.tex > mainSupplement_todo.tex
	pdflatex mainSupplement_todo.tex
	bibtex mainSupplement_todo.aux
	pdflatex mainSupplement_todo.tex
	pdflatex mainSupplement_todo.tex

mainSupplement.bbl: mainSupplement.aux
	bibtex mainSupplement.aux

mainSupplement.aux: mainSupplement.tex
	pdflatex mainSupplement.tex


otherFigures.pdf: otherFigures.tex ${supfigures} ${supfigurespdf} ${suptex} supplementReset.tex
	pdflatex otherFigures.tex
	pdflatex otherFigures.tex

pdfFigures:
	cd figures/; pdflatex Fig1.tex
	cd figures/; pdflatex Fig2.tex
	cd figures/; pdflatex Fig3.tex
	cd figures/; pdflatex Fig4.tex
	cd figures/; pdflatex Fig5.tex
	cd figures/; pdflatex Fig7.tex
	cd figures/; pdflatex Fig2Africa.tex
	cd figures/; pdflatex Fig2SupMix3.tex
	cd figures/; pdflatex Fig2Asia.tex
	cd figures/; pdflatex Fig2k3.tex
	cd figures/; pdflatex Fig2k4.tex
	cd figures/; pdflatex Fig3Sup.tex
	cd figures/; pdflatex Fig2inVitro.tex

tmpmain.pdf: tmpmain.tex
	pdflatex tmpmain.tex

tmpmain.tex: main.tex Makefile
	sed -e "/Manuscript submitted to eLife/d" elife.cls > tmpelife.cls
	sed -e "s/\\\documentclass\[9pt,lineno\]{elife}/\\\documentclass\[9pt\]{tmpelife}/" main.tex > tmpmain.tex

plain.pdf: plain.tex tmpmain.pdf
	pdflatex plain.tex
	bibtex plain.aux
	pdflatex plain.tex
	pdflatex plain.tex

plain.tex: main.tex Makefile
	sed -e "s/elife/article/" \
	 -e "s/\\\usepackage{hyperref}/\\\usepackage{hyperref, natbib, fullpage}\n/" \
	 -e "/\\\maketitle/d" \
	 -e "s/\\\contrib\[$$\\\dagger$$\]{These authors contributed equally to this work}/\\\maketitle\n\\\footnotetext[2]{These authors contributed equally to this work}/" \
	 -e "s/\\\corr{gil.mcvean@bdi.ox.ac.uk}{GM}/\\\footnotetext[1]{For correspondence: gil.mcvean@bdi.ox.ac.uk}/" \
	 -e "s/width=\\\textwidth/width=0.8\\\textwidth/g" \
	 -e "/figsupp/d" \
	 -e "s/tabledata/caption/" \
	 -e "s/\\\end{document}/\\\includepdf[pages=21-27]{tmpmain.pdf}\n\\\end{document}/" \
	 -e "s/\\\bibliography{mixedIBD.bib}/\\\bibliographystyle{chicagoa}\n\\\bibliography{mixedIBD.bib}/" \
	  main.tex > plain.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out main.pdf

