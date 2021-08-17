# A prototype Makefile for H1
#
# last update: $Date: 2010/10/08 13:39:31 $
# by         : $Author: msteder $
#
# The name of this module
MODULENAME = H1Examples
#
# macro definitions
# DEFINES = 
#
# debug on or off (comment next line to disable debugger information)
DEBUG = -g -Wall
#
# The shared library to make (if any): e.g. libH1Tracks.so
# SLIB_TARGETS = 
SLIB = libAzimuthalAnalyzer.so
#
# The static library to make (if any): e.g. libH1Tracks.a
# LIB_TARGETS = 
LIB = libAzimuthalAnalyzer.a
#
# The binaries you would like to make, you have to provide the final linking
BIN = AzimuthalAnalyzer

#ifeq ($(STANDALONE),yes)
  # no rerun_rec
#else
#  BIN += rerun_rec 
#endif

#
# Objects you need for the libraries
#OBJS = AzimuthalAnalyzer.o
OBJS = FidVolCut.o

# OBJS += UserEvent.o UserFiller.o

# Dictionary objects, default is: DICTOBJS = $(OBJS:.o=Dict.o)
DICTOBJS = $(OBJS:.o=Dict.o) 
DICTOBJS += TDetectQedcDict.o

# Now your special libs and includes
# EXTRA_INCLUDEDIRS = ~/my_favourite_library/include
# EXTRA_LIBDIRS = ~/my_favourite_library/lib
# EXTRA_LINKFLAGS = -lmy_favourite_library
# Dummy objects
# DUMMYOBJS =

# Dummy dictionary objects
# DUMMYDICTOBJS = $(DUMMYOBJS:.o=Dict.o)


#dict: rootcint -f AzimuthalAnalyzer.C -c $(INCLUDES)


#extra files to be removed when calling 'make clean' or 'make clobber'
#EXTRA_CLEAN_RULES=
#EXTRA_DISTCLEAN_RULES=

# defines LIBS
# defines Rules for .c .C .cpp .f Dict.C .so
ifeq ($(H1BASICMK),local)
  include  ../oo_tools/H1Basic.mk
else
  include  $(H1DIST)/releases/$(H1CURRENT)/oo_tools/H1Basic.mk
endif                                                             

#ROOTLIBS=$(shell root-cofig --libs)

# include the dependencies if created via 'make depend'
-include Makefile.dep

# the executables
AzimuthalAnalyzer: $(BUILDDIR)/AzimuthalAnalyzer.o $(LIB) outputdirs
  @echo ""
  @echo Linking test program for ods access: $@ ...
  $(LOADER) \
  $(BUILDDIR)/AzimuthalAnalyzer.o \
  -lAzimuthalAnalyzer \
  $(OO_LIBS) \
  -lH1UserCim -lH1Clusters -lH1SubDetInfo \
  -lH1Tools -lH1Arrays -lH1Hat -lH1Mods \
  -lH1Ods -lH1MagfieldOO \
  -lH1Geom -lH1Skeleton \
  -lH1HadronicCalibration \
  -lH1OOBanks -lH1Pointers -lH1Skeleton  \
  -lH1Steering \
  -lH1Tracks -lH1PhysUtils  -lH1JetFinder -lH1Finder -lH1Clusters \
  -lH1Geom -lH1Ods -lH1OOBanks -lH1Tracks -lH1Clusters -lH1Mods -lH1Hat -lH1Skeleton -lH1Steering \
  -lH1Pointers -lH1Arrays -lH1Tools -lH1Benchmarks -lH1Arrays -lH1Binning \
  -lH1PhysUtils  -lH1JetFinder -lH1Finder -lH1HadronicCalibration -lH1CaloTrigger \
  -lH1Calculator -lH1CalcPointers -lH1CalcWeights -lH1SoftLeptonId  -lH1Analysis \
  -lH1MagfieldOO -lH1Cuts -lH1TrkFinder -lH1Selection -lH1ElecCalibration -lH1QCDFunc  -l${LHAPDFLIB} \
  \
  -L$(H1_LIBS) -L$(CERN_LIBS) \
  -L../lib/x86_64-centos7-gcc9-opt \
  $(ROOTLIBS) -lEG -l${FORTRANLIB} \
  -Wl,-rpath=../lib/amd64_rhel60 \
  -o $(BINDIR)/AzimuthalAnalyzer
  @echo ""
  @echo "copying .pcm files to libdir"
  cp *.pcm $(LIBDIR)
  @echo ""
  @echo $@ created ...
  @echo ""