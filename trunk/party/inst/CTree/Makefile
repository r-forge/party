DOCUMENT=conditionaltrees

TEXT= \
	$(DOCUMENT).tex \
	head.tex \
	abstract.tex \
	introduction.tex \
	methods.tex \
	applications.tex \
	equivalence.tex \
	discussion.tex \
	appendix.tex \
	benchdata.tex \
	acknowledgements.tex \
	$(DONE)

PICS= \
	applications-gbsg2.pdf \
	applications-glaucoma.pdf \
	applications-mammo.pdf \
	nmiratiodens.pdf \
	nmisim.pdf \
	power_split.pdf \
	ratiobox.pdf \
	ratioboxQUEST.pdf \
	ratioboxQUEST_GUIDE.pdf \
	$(DONE)

BIBS = \
	$(DOCUMENT).bib \
	$(DONE)

all	: $(DOCUMENT).pdf


applications.tex:	applications.Rnw 
	echo "library(utils); Sweave(\"applications.Rnw\");" | R --vanilla 

equivalence.tex:	equivalence.Rnw 
	echo "library(utils); Sweave(\"equivalence.Rnw\");" | R --vanilla 


thpdf	: $(DOCUMENT).pdf
	thumbpdf $(DOCUMENT).pdf
	R CMD texi2dvi --pdf $(DOCUMENT).tex

$(DOCUMENT).pdf	: $(TEXT) $(BIBS)
	R CMD texi2dvi --pdf $(DOCUMENT).tex

distclean   :
	-rm -f `ls  $(DOCUMENT).* | sed -e /\.tex/d`
	-rm -f *~

clean   :
	-rm -f `ls  $(DOCUMENT).* | sed -e /\.tex/d`
	-rm -f *~

