CC        = gcc
CFLAGS    = -Wall -g -O3 -Iwrappers/
LFLAGS    = -Wall -g -O3 
OBJS      = cupidconfigD.o cupidclumpfind.o test.o 
EXE       = test wrappers/mers.o wrappers/ast.o

all : $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o $(EXE)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

