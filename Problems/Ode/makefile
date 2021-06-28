CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt cos.png epidemic.png threeBody.png


out.txt: main
	./$< > $@

main: main.o gravity.o
gravity.o: gravity.c

cos.png: cos.txt
	echo '\
	set term png;\
	set output "$@";\
	set key bottom right;\
	set fontsize 0.7;\
	set xlabel "t";\
	set ylabel "y";\
	set title "ODE of cos(t)";\
	plot \
	"$<" using 1:2 with lines color red lw 3 title "y from ODE" \
	,"$<" using 1:3 with lines color blue lw 3 title "dy/dt from ODE" \
	' | pyxplot

epidemic.png: epidemic.gpi Epidemic2.txt Epidemic5.txt Epidemic10.txt
	pyxplot $<

threeBody.png: threeBody.gpi threeBodyProblem.txt
	pyxplot $<


.PHONEY:clean
clean:
	$(RM) *txt main *png
