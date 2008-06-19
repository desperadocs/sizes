#   UNIX Makefile for the program: sizes
# author:  Pete R. Jemian
#          Late-Nite(tm) Software
#          jemian@anl.gov
# revised: 12 August 1996
# revised: 3 December 2004

SOURCES  = sizes.c 
SOURCES += grid.c
SOURCES += lu.c
SOURCES += maxent.c
SOURCES += recipes.c
SOURCES += regular.c
SOURCES += stats.c
SOURCES += nnls2.c

INCLUDES = grid.h recipes.h nnls2.h

OBJECTS = $(SOURCES:%.c=%.o)

# for use with debugger
CFLAGS += -g

# MALLOC_CHECK_ = 1

ARCHIVE_FILES  = Makefile
ARCHIVE_FILES += $(SOURCES) 
ARCHIVE_FILES += $(INCLUDES) 
ARCHIVE_FILES += grid.h recipes.h nnls2.h
ARCHIVE_FILES += sizes sizesgui.tcl
ARCHIVE_FILES += test.dis test.sas
ARCHIVE_FILES += test-fm.cmd
ARCHIVE_FILES += test-fr.cmd
ARCHIVE_FILES += docs

# no _real_ preference for GNU C compiler here.  You can change it.
CC = gcc


# You're going to need the math library.
LD_LIBS += -lm
# use Electric Fence to check memory errors (ONLY) when debugging
#LD_LIBS += -lefence

sizes: $(OBJECTS) Makefile
	$(CC) -o sizes $(OBJECTS) $(LD_LIBS)

all: sizes archive

test: sizes
	sizes test-f.cmd

all-tests: sizes
	sizes test-f.cmd
	sizes test-i.cmd
	sizes test-n.cmd
	sizes model2.cmd
	sizes ps-dsm.cmd
	sizes ps-smr.cmd

clean: clean_core
	/bin/rm -f *.o *%
	@cd ./docs && make realclean

clean_core:
	/bin/rm -f core core.*

archive tar: sizes.tar.gz

guide :: 
	@cd ./docs && make

sizes.tar.gz: $(SOURCES) $(ARCHIVE_FILES) guide
	rm -f sizes.tar.gz
	tar cf sizes.tar $(ARCHIVE_FILES)
	gzip sizes.tar
