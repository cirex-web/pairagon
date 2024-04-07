###############################################################################
# Makefile for zoe
###############################################################################

###########
# LIBRARY #
###########

LIBRARY = lib/libzoe.a

OBJECTS = \
	src/zTools.o\
	src/zTBTree.o\
	src/bntree/bntree.o\
	src/bntree/cache.o\
	src/bntree/celllistitem.o\
	src/bntree/intlistitem.o\
	src/bntree/node.o\
	src/zAlignment.o\
	src/zAlnFeature.o\
	src/zConseq.o\
	src/zDistribution.o\
	src/zDNA.o\
	src/zDuration.o\
	src/zEstseq.o\
	src/zEstseqSearch.o\
	src/zFastaFile.o\
	src/zFeatureFactory.o\
	src/zGTF.o\
	src/zHardCoding.o\
	src/zHMM.o\
	src/zHMM_State.o\
	src/zMath.o\
	src/zMathTables.o\
	src/zModel.o\
	src/zPhasePref.o\
	src/zProtein.o\
	src/zScanner.o\
	src/zSeedUtils.o\
	src/zSequence.o\
	src/zSfeature.o\
	src/zStopSeq.o\
	src/zTransition.o\
	src/zTrellis.o\
	src/zViterbi.o\
	src/zPairTransition.o\
	src/zPairTrellis.o\
	src/zPairViterbi.o\

PERL_LIBS = \
	lib/perl5/Alignment.pm

PAIRAGON_SCRIPTS= \
	bin/Pairagon.pl \
	bin/alignmentConvert.pl

#################
# TEST PROGRAMS #
#################

LFLAGS = -L. -Llib -lzoe -lm 

EXE15 = bin/pairagon
SRC15 = src/pairagon.c
OBJ15 = $(SRC15:.c=.o)

EXE16 = bin/pairameter_estimate
SRC16 = src/pairameter_estimate.c
OBJ16 = $(SRC16:.c=.o)

EXECUTABLES = $(EXE15) $(EXE16)

DESTDIR ?= /usr

TABLE = src/zMathTables.c
TABLEGEN = src/zFloatwiseScoreAdd

########################
# Program Options      #
########################

BUILD = 20090501
PAIRAGON_VERSION = 1.1
PAIRAGON_PROGRAM = Pairagon
GCC_FLAGS = "-W -Wall -Werror -ansi -pedantic -Wwrite-strings"
PAIRAGON_DIST_FLAGS = "-DBUILD=\\\"$(PAIRAGON_PROGRAM)\ version\ $(PAIRAGON_VERSION)\ build\ $(BUILD)\\\" -DNDEBUG -O2"

# In order to make the include file to contain header file dependencies,
# gcc needs a flag to construct this.  In another funny SunOS stunt, this
# is entirely different in their cc.
INCLUDE_FLAGS = "-MM"

SPARC_INCLUDE_FLAGS= "-xM1"

ALPHA_FLAGS = "-c99 -msg_enable level5"

SPARC_FLAGS = "-O "

MACOSX_FLAGS = "-Wall -W -O2 -g"
MACOSX_G5_FLAGS="-fast"
# You might have to use -mcpu=7450 if it is not a G5 processor!

############################
# End Program Options      #
############################

################################################################################
# TARGETS 
################################################################################


all: include $(LIBRARY) $(EXECUTABLES) 

$(LIBRARY): $(OBJECTS)
	rm -f $(LIBRARY)
	ar cr $(LIBRARY) $(OBJECTS)
	ranlib $(LIBRARY)

include:
	$(CC) $(INCLUDE_FLAGS) src/*.c > Makefile.include
#	src/get-glib-flags.sh >> Makefile.include

pairagon-exe: $(LIBRARY) $(EXECUTABLES)

clean:
	rm -f src/*.o src/bntree/*.o $(LIBRARY) $(EXECUTABLES) $(TABLE) $(TABLEGEN)

product:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=i686 CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

gcc-distribution:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=i686 CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

linux-generic:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

linux:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=i686 CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

linux64:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=nocona CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

alpha:
	make all CC="cc" CFLAGS+=$(ALPHA_FLAGS) CFLAGS+=$(DIST_FLAGS)

macosx:
	make all CC="gcc" CFLAGS+=$(MACOSX_FLAGS) CFLAGS+=$(DIST_FLAGS)

macosxg5:
	make all CC="gcc" CFLAGS+=$(MACOSX_FLAGS) CFLAGS+=$(MACOSX_G5_FLAGS) CFLAGS+=$(DIST_FLAGS)

sparc:
	make all CC="cc" INCLUDE_FLAGS=$(SPARC_INCLUDE_FLAGS) CFLAGS+=$(SPARC_FLAGS) CFLAGS+=$(DIST_FLAGS)

pairagon-linux:
	make pairagon-exe CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=$(PAIRAGON_DIST_FLAGS) $(GLIB_FLAGS)

pairagon-sparc:
	make pairagon-exe CC="cc" CFLAGS+=$(SPARC_FLAGS) CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairagon-macosx:
	make pairagon-exe CC="cc" CFLAGS+=$(MACOSX_FLAGS) CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairagon-distro:
	make pairagon-exe CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=i686 CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairameter-estimate: $(EXE16)

install: 
	mkdir -p $(DESTDIR)/bin
	install -m 0755 $(EXECUTABLES) $(DESTDIR)/bin
	for f in $(PAIRAGON_SCRIPTS); do \
		install -m 0755 $$f $(DESTDIR)/$$f ;\
	done
	mkdir -p $(DESTDIR)/lib/perl5
	install -m 0644 $(PERL_LIBS) $(DESTDIR)/lib/perl5

#######################
# Excecutable Section #
#######################

$(EXE15): $(OBJ15) $(LIBRARY)
	$(CC) -o $(EXE15) $(CFLAGS) $(OBJ15) $(LFLAGS) $(GLIB_LFLAGS)

$(EXE16): $(OBJ16) $(LIBRARY)
	$(CC) -o $(EXE16) $(CFLAGS) $(OBJ16) $(LFLAGS) $(GLIB_LFLAGS)
###################
# Inference Rules #
###################

.SUFFIXES : .c .o

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<


#######################
# Source Dependencies #
#######################

$(TABLEGEN): src/zFloatwiseScoreAdd.c
	$(CC) $(CFLAGS) -o $(TABLEGEN) src/zFloatwiseScoreAdd.c -lm

$(TABLE): $(TABLEGEN)
	./$(TABLEGEN) > $(TABLE)

include Makefile.include

