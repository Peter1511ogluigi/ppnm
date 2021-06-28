
CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
IMGFLS  =  plot.png

.PHONEY: default
default: out.txt data.txt $(IMGFLS)

#plot.png: plot.gpi data.txt
#	pyxplot $<


plot.png: plot.gpi data.txt
		pyxplot $< 1> /dev/null
plot.gpi: data.txt makefile
		@echo 'set terminal png dpi 600'						>  $@
		@echo 'set output "plot.png"'						>> $@
		@echo 'set key outside' 								>> $@
		@echo 'set yrange [-1.5:1.5]'							>> $@
		@echo 'set xlabel "$$x$$"'								>> $@
		@echo 'set ylabel "$$\psi$$"'								>> $@
		@echo 'set title "Quantum particle in a box (with dimension 100x100)"'	>> $@
		@echo 'plot \\'																>> $@
		@echo ' "data.txt" using 1:2 w l lw 6 color black  title "$$n=1$$ (analytical)" 		\\'>> $@
		@echo ',"data.txt" using 1:3 w l lw 3 color red    title "$$n=1$$ (calculated)" 	\\'>> $@
		@echo ',"data.txt" using 1:4 w l lw 6 color black  title "$$n=2$$ (analytical)" 		\\'>> $@
		@echo ',"data.txt" using 1:5 w l lw 3 color cyan  title "$$n=2$$ (calculated)" 	\\'>> $@
		@echo ',"data.txt" using 1:6 w l lw 6 color black  title "$$n=3$$ (analytical)" 		\\'>> $@
		@echo ',"data.txt" using 1:7 w l lw 3 color limegreen title "$$n=3$$ (calculated)" 	\\'>> $@



data.txt: main
		./$< > $@

out.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) *.o *.txt $(IMGFLS)
