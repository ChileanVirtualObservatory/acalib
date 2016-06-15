CC        = gcc
CFLAGS    = -Wall -g -O3 -Iwrappers/
LFLAGS    = -Wall -g -O3 
OBJS      = cupidconfigI.o cupidconfigD.o cupidclumpfind.o test.o 
EXE       = test wrappers/mers.o wrappers/ast.o
SUBDIRS   = wrappers


all : $(OBJS) $(SUBDIRS)
	$(CC) $(LFLAGS) $(OBJS) -o $(EXE)


subdirs: $(SUBDIRS)

.PHONY: subdirs $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@


%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

