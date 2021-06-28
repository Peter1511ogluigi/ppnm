CFLAGS= -Wall -O -std=gnu11
CFLAGS+=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs)

.PHONEY: default
default: out.txt plots.png



out.txt: main
	./$< > $@



plots.png: plots.pgi points.txt
		pyxplot $<

##functions.txt: main
##	./$< > $@

main: main.o ann.o qnewton.o
ann.c: qnewton.c
.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
