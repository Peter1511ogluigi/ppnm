CFLAGS= -Wall -O -std=gnu11
CFLAGS+=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs)

.PHONEY: default



out.txt: main
	./$< > $@

##functions.txt: main
##	./$< > $@

main: main.o adaptiveIntegraton.o 
.PHONEY: clean
clean:
	$(RM) main *.o *.txt
