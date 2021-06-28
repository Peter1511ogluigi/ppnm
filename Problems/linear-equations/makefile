CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: out.txt timeplot.png


out.txt: main
	./main > out.txt

main.o: main.c


timeplot.png: timeplot.gpi timeData.txt
		pyxplot $<

main: main.o diffClock.o
main.c: diffClock.c


.PHONEY: clean


clean:
	$(RM) *.o out* main *.png
