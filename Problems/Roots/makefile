CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: out.txt Hydrogen.png


out.txt: main
	./$< > $@

Hydrogen.txt: main
		./$< > $@

main: main.o roots.o odesolver.o


Hydrogen.png: Hydrogen.pgi Hydrogen.txt
	pyxplot $<




.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
