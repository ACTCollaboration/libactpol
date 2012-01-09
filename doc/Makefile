
DOC = pointing

.PHONY: ps pdf all clean default

default: pdf

ps:
	latex $(DOC)
	bibtex $(DOC)
	latex $(DOC)
	latex $(DOC)
	dvips $(DOC) -o

pdf:
	latex $(DOC)
	bibtex $(DOC)
	latex $(DOC)
	latex $(DOC)
	dvipdf $(DOC)

all: ps pdf

clean:
	-rm *.aux *.log *.dvi *.bbl *.blg

