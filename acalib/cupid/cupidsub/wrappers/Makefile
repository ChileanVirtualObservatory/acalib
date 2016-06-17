CC        = gcc
PYINCLUDE = /usr/include/python3.5m/
CFLAGS    = -Wall -g -O3 -Iincludes -I$(PYINCLUDE)
OBJS      = hds.o mers.o ast.o ndf.o

all : $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

