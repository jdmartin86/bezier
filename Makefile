#
# make file for bezier 
#

CFLAGS = -g -std=c99 -O -fopenmp -o
INCLUDE = -lm 

# rules
all: bezier.c
	gcc $(CFLAGS) bezier $< $(INCLUDE)

clean:
	rm -f bezier Makefile~ *.txt *.txt~ *.py~ *.c~ *.h~  
