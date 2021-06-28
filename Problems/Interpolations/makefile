
CFLAGS 		 +=  $(shell gsl-config --cflags)
LDLIBS 		 +=  $(shell gsl-config --libs)
LDLIBS 		 +=  -lm -Wall

CC		  		=  gcc
PIPE     		=  cat

INPUT   		=  cos_tabulated.txt
LINOUT   		=  lin_spline_output.txt
QUADOUT  		=  quad_spline_output.txt
CUBICOUT  	=  cubic_spline_output.txt
EXE   	  	=  main
OBJ  	  		=  main.o linSpline.o integrateFunction.o quadSpline.o cubicSpline.o utilities.o
HDRS    	  =  linSpline.h integrateFunction.h quadSpline.h cubicSpline.h utilities.h
IMGFLS  	  =  linear.png quad.png cubic.png
TXT     		=  out.txt

.PHONEY: default
default: $(IMGFLS) print

print: $(TXT)


linear.png: linear.gpi lin_spline_output.txt cos_tabulated.txt
	pyxplot $<

quad.png: quad.gpi quad_spline_output.txt cos_tabulated.txt
	pyxplot $<

cubic.png: cubic.gpi cubic_spline_output.txt cos_tabulated.txt
	pyxplot $<

$(LINOUT) $(QUADOUT) $(CUBICOUT) $(TXT): $(EXE) $(INPUT)
	./$< $(INPUT) $(LINOUT) $(QUADOUT) $(CUBICOUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(OUTPUT) $(OBJ) $(IMGFLS) $(TXT) $(LINOUT) $(QUADOUT) $(CUBICOUT)
