#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

erf.gnuplot.png: test.gpi
	cat test.gpi | sed 's/pyxplot/gnuplot/' | gnuplot

erf.pyxplot.png: data.txt makefile tabulated.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

gamma.gnuplot.png: gamma.gpi
		cat gamma.gpi | sed 's/pyxplot/gnuplot/' | gnuplot

gamma.pyxplot.png: data.txt makefile tabulated.data.txt
		echo '\
		set terminal png;\
		set output "$@";\
		set key top left;\
		set tics out;\
		set xlabel "x";\
		set ylabel "y";\
		set title "Gamma function";\
		plot \
		 "$<" using 1:2 with line title "gamma from math.h"\
		,"$<" using 1:3 with line title "gsl-sf-gamma"\
		,"$<" using 1:4 with line title "homemade Gamma"\
		,"tabulated.data.txt" using 1:2 with points title "tabulated data"\
		' | tee gamma.gpi | pyxplot



test.gpi: erf.pyxplot.png
gamma.gpi: gamma.pyxplot.png


data.txt: main
	./$< > $@

main: main.o erf.o


.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
