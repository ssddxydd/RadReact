.SUFFIXES:
.SUFFIXES: .out .o .c
# Set the path and parameters below to the correct values when needed
INCDIR = /usr/local/include
LIBDIR = /usr/local/lib
CC=gcc
CFLAGS= -I$(INCDIR) -L$(LIBDIR)

NRLIBS=-lm
#.c.o: ; $(CC) $(CFLAGS) -c $<
#.o.out: ; $(CC) $(CFLAGS) -g -o $@ $^ $(NRLIBS) 

.c.o: ; $(CC) $(CFLAGS) -O3 -c $<
.o.out: ; $(CC) $(CFLAGS) -O3 -o $@ $^ $(NRLIBS) 

all : main.out

main.out : main.o  bsstep.o \
	rkqs.o rkck.o \
	nrutil.o odeint.o\
	mmid.o pzextr.c\
	beschb.o bessik.o chebev.o\
	fairy.o qromo.o polint.o\
	midpnt.o simpr.o\
	ludcmp.o

main.o : main.c
beschb.o : beschb.c
bessik.o : bessik.c
bsstep.o : bsstep.c
chebev.o : chebev.c 
rkqs.o : rkqs.c 
rkck.o : rkck.c 
nrutil.o : odeint.c 
odeint.o : odeint.c 
mmid.o : mmid.c 
pzextr.o : pzextr.c
polint.o : polint.c
fairy.o : fairy.c
qromo.o : qromo.c
midpnt.o : midpnt.c
simpr.o : simpr.c
ludcmp.o : ludcmp.c

clean : ;
	/bin/rm *.o *.out
