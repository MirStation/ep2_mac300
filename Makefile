# Makefile
FLAGS = -g -Wall -std=gnu99 -pedantic -lm -O3
TARGETS = cgm sspdsgen
OBJECTS = cgm.o sspdsgen.o

all : ${TARGETS}

cgm : cgm.o
	gcc -o $@ $? ${FLAGS}

sspdsgen : sspdsgen.o
	gcc -o $@ $? ${FLAGS}

%.o : %.c
	gcc -c -o $@ $? ${FLAGS} -I.

clean :
	rm -f *~ ${OBJECTS} ${TARGETS}

tar:
	tar cvfz ep2.tgz Makefile *.c *.h *.pdf
