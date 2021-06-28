CC = gcc
CFLAGS = -std=gnu11 -O -Wall -fwrapv
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: main
	./$< > $@


main: main.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	$(RM) *.o out* main
