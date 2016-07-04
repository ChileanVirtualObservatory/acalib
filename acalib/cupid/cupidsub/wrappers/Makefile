CC        = gcc
PYINCLUDE = /home/martin/miniconda2/include/python2.7
CFLAGS    = -Wall -g -fPIC -O3 -Iincludes -I$(PYINCLUDE)
#OBJS      = hds.o mers.o ast.o ndf.o
OBJS      = mers.o ast.o 

all : $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

