CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt fit.png


fit.png: fit.pgi exp_vals.txt linfit.txt
	pyxplot $<

out.txt: main
	./main > out.txt


main: main.o leastSquare.o gramSchmidt.o

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main.c: leastSquare.c

leastSquare.c: gramSchmidt.c


clean:
	$(RM) *.o out* main
