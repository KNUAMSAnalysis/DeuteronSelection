CXX = g++
CXXFLAGS = -W -Wall -Wno-unused-parameter -Wno-unknown-pragmas -DAMS_ACQT_INTERFACE -D_PGTRACK_ -Wno-unused-variable -Wno-extra

# AMS Global Environment
CVMFS_AMS_OFFLINE = /cvmfs/ams.cern.ch/Offline

############################# CERN libraries ##################################
CERNLIBS = -lgeant321 -lpacklib -lmathlib -lkernlib
CERN_LEVEL = 2005.gcc64

ifndef CERNDIR
CERNDIR = $(CVMFS_AMS_OFFLINE)/CERN/$(CERN_LEVEL)
endif

CERNSRCDIR = $(CERNDIR)

ifndef AMSLIB
AMSLIB = /$(CVMFS_AMS_OFFLINE)/lib/linux/gcc64
endif

ifndef NAGDIR
NAGDIR = $(CVMFS_AMS_OFFLINE)/CERN/NagLib
endif
######################### End of CERN library settings ########################

# AMS Offline Software Related Includes
INCLUDES = -I${ROOTSYS}/include -I${AMSWD}/include -I./include -I/afs/cern.ch/user/w/wyjang/AMS-ACsoft/io
NTUPLE_PG = $(AMSWD)/lib/linuxx8664gcc5.34/ntuple_slc6_PG.so
############ End of AMS Offline Software related includes

# ROOT Related Settings
ROOTLIBS = $(shell root-config --libs) -lASImage -lRIO -lNet -lNetx -lMinuit -lTMVA -lMLP -lXMLIO -lTreePlayer
############ End of ROOT related settings

# ACSOFT Flags
#ACSOFTFLAGS = `acsoft-config --definitions` `acsoft-config --cflags` `acsoft-config --auxcflags` `acsoft-config --libs` `acsoft-config --auxlibs`
ACSOFTFLAGS = `acsoft-config --definitions` `acsoft-config --cflags` `acsoft-config --cflags-include` `acsoft-config --auxcflags-include` `acsoft-config --libs` `acsoft-config --auxlibs`

TARGET = bin/main

all : $(TARGET)

$(TARGET) : obj/main.o obj/selector.o
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -o $@ $^ $(NTUPLE_PG)

obj/selector.o : src/selector.cxx
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c -o $@ $^

obj/main.o : src/main.cxx
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) -Wno-extra $(ROOTLIBS) $(INCLUDES) -c -o $@ $^

clean :
	rm -rf obj/*.o bin/main

