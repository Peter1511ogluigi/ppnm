CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt timePlot.png


out.txt: main
	./main > out.txt

timePlot.png: timePlot.gpi timeData.txt
			pyxplot $<

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

jacobian.o: jacobian.c


main: main.o jacobian.o



.PHONEY: clean


clean:
	$(RM) *.o out* main *.txt *.png
