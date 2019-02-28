hd= /home/tinker


FC= ifc
CC= gcc
FLINK= ifc
CLINK= gcc 
CFLAGS=  $(FFTW_CFLAGS)
FLIB= -L$(hd)/cosmo/lib -lfutil -lcutil 
CLIB= -L${HOME}/cosmo/lib -lcutil -lm 
XDIR= $(HOME)/exec


#CC = gcc
#CLINK = gcc

OBJS1=	lh_optimize.o amoeba.o amotry.o powell.o linmin.o brent.o mnbrak.o f1dim.o
lh_optimize:	$(OBJS1)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)
	cp $@ $(HOME)/exec/

OBJS2=	lh_mcmc.o gasdev.o ran1.o
lh_mcmc:	$(OBJS2)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS2) $(CLIB)
	cp $@ $(HOME)/exec/

OBJS3= hypercube_eigenspace.o jacobi.o gaussj.o
hypercube_eigenspace:	$(OBJS3)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS3) $(CLIB)
	cp $@ $(HOME)/exec/

clean:
	rm -f *.o
