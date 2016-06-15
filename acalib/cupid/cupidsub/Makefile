CC        = gcc
CFLAGS    = -Wall -g -O3 -Iwrappers/includes
LFLAGS    = -Wall -g -O3 -lm 
OBJS      = cupidndfclump.o cupidcfaddpixel.o cupidcfclump.o cupidcfdeleteps.o cupidcferode.o cupidcfidl.o cupidcfnebs.o cupidcfscan.o  cupidcfxfer.o cupidcfmakeps.o cupidcffreeps.o cupiddefminpix.o cupidconfigrms.o cupidcflevels.o cupidconfigI.o cupidconfigD.o cupidclumpfind.o test.o  wrappers/mers.o wrappers/ast.o wrappers/hds.o
EXE       = test 
SUBDIRS   = wrappers


all : $(OBJS) $(SUBDIRS)
	$(CC) $(OBJS) $(LFLAGS)-o $(EXE)


subdirs: $(SUBDIRS)

.PHONY: subdirs $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@


%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXE)

