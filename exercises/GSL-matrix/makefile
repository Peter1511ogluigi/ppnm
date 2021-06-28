
#CFLAGS = -Wall -O1 -std=gnu11
#CFLAGS += $(shell /usr/bin/gsl-config --cflags)
#LDLIBS += $(shell /usr/bin/gsl-config --libs)

CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
CXXFLAGS = -Wall -std=c++0x `gsl-config --cflags`
LDLIBS = -lstdc++ `gsl-config --libs`

out.txt: main
	./$< > out.txt

.PHONEY: clean
clean:
	$(RM) main.o main out.txt
