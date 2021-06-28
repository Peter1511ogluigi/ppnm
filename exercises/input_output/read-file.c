#include <stdio.h>
#include <math.h>

int main(){
int items; 
double x;
FILE* my_out_stream = fopen("input.txt","r");					 //input.txt is the txt file you want to have knowledge of (want to know sine/cosine)
while(items!=EOF){							// making a loop which reads each line from the file, and stops when the last line has been read.
		items = fscanf(my_out_stream,"%lg",&x);
		printf("x=%g sin(x)=%g cos(x)=%g\n",x,sin(x),cos(x));	
	}
	
fclose(my_out_stream);
return 0;
}
