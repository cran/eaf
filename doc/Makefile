
R_PROG = R
R_HOME = `${R_PROG} RHOME`

.PHONY: pdf cleanall clean

pdf: info.tex
	pdflatex info.tex
	pdflatex info.tex
	$(MAKE) clean

%.tex: %.Rnw
	echo "library(tools); Sweave('$<')" | ${R_PROG} --no-save --no-restore --no-site-file --no-init-file --silent

%.R: %.Rnw
	echo "library(tools); Stangle('$<')" | ${R_PROG} --no-save --no-restore --no-site-file --no-init-file --silent


cleanall:
	@$(RM) info.tex *.log *.aux *.dvi *.ps *.pdf figs/* *~ */*~

clean:
	@$(RM) info.tex *.log *.aux *.dvi *.ps *~  Rplots.pdf fig-*.pdf

