CC        = gcc
PYINCLUDE = /home/martin/miniconda2/include/python2.7
CFLAGS    = -Wall -g -fPIC -O3 -Iwrappers/includes -I$(PYINCLUDE)
LFLAGS    = -Wall -g -O3 -lm -lpython2.7 
#OBJS      = cupidndfclump.o cupidcfaddpixel.o cupidcfclump.o cupidcfdeleteps.o cupidcferode.o cupidcfidl.o cupidcfnebs.o cupidcfscan.o  cupidcfxfer.o cupidcfmakeps.o cupidcffreeps.o cupiddefminpix.o cupidconfigrms.o cupidcflevels.o cupidconfigI.o cupidconfigD.o cupidclumpfind.o test.o  wrappers/mers.o wrappers/ast.o wrappers/hds.o wrappers/ndf.o
OBJS      = cupidcfaddpixel.o cupidcfclump.o cupidcfdeleteps.o cupidcferode.o cupidcfidl.o cupidcfnebs.o cupidcfscan.o  cupidcfxfer.o cupidcfmakeps.o cupidcffreeps.o cupiddefminpix.o cupidconfigrms.o cupidcflevels.o cupidconfigI.o cupidconfigD.o pycupidclumpfind.o test.o  wrappers/mers.o wrappers/ast.o 
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

