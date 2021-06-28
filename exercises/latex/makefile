CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

latex.pdf: latex.tex fig-pyxplot.pdf fig-pyxplot-2.pdf
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

out.data.txt:main
	./$< > $@

main:main.o myexponential.o

fig-pyxplot.pdf: fig-pyxplot.gpi out.data.txt
	pyxplot $<
fig-pyxplot-2.pdf: fig-pyxplot-2.gpi out.data.txt
	pyxplot $<

clean:
	$(RM) latex.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex
