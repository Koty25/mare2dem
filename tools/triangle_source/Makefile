

#-------------------------------------------------------------------------------
# Specify your C Compiler
#-------------------------------------------------------------------------------

# Gcc compiler:
#CC = cc
#COPTS   = -O2

# Intel Compiler:
CC = icc
COPTS   = -O2  -fp-model precise -fp-model source # -fast is best optimization on Mac OS X

#-------------------------------------------------------------------------------
# Compilation Commands
#-------------------------------------------------------------------------------

TARGETS		=  clean triangle 

all:		$(TARGETS)
 
clean: 	
	rm -f *.o triangle
 
triangle: triangle.c
		$(CC) $(COPTS)  -o $@ $<
		
 