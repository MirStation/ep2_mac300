# Makefile
FLAGS = -g -Wall -ansi -pedantic -lm -O3
TARGETS = cgm 
OBJECTS = cgm.o

all : ${TARGETS}

cgm : cgm.o
	gcc -o $@ $? ${FLAGS}

%.o : %.c
	gcc -c -o $@ $? ${FLAGS} -I.

clean :
	rm -f *~ ${OBJECTS} ${TARGETS}

tar:
	tar cvfz ep2.tgz Makefile *.c *.h *.pdf
