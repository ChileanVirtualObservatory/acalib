CC        = gcc
CFLAGS    = -Wall -g -O3 -Iincludes
OBJS      = hds.o mers.o ast.o ndf.o

all : $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

