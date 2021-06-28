CFLAGS= -Wall -O -std=gnu11
CFLAGS+=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs)

.PHONEY: default
default: out.txt plot.png



out.txt: main
	./$< > $@

plot.png: plot.gpi data.txt
		pyxplot $<


##functions.txt: main
##	./$< > $@

main: monteCarlo.c
.PHONEY: clean
clean:
	$(RM) main *.o *.txt
