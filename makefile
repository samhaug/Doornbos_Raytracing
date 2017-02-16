# Make file for programs in this directory.  Object location is assumed to
#   be in the superior directory to this.
#                 
DEBUG=	-O    
FC= /usr/bin/gfortran
FFLAGS= ${DEBUG}
CFLAGS=	-g -O

PROGS=  ray1d \

ray1d: ray1d.o rayps.o modps.o 
	${FC} ${FFLAGS} -o ./ray1d ray1d.o rayps.o modps.o

clean:

	/bin/rm -f *.o ${PROGS}
