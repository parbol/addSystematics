ObjSuf        = o
SrcSuf        = cc
ExeSuf        = cpp
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit
ROOTLIBS      = -lGenVector -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTMVA
#ROOTGLIBS     = -lGenVector -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES) -I ${LHAPDFINC}


LD		= g++
LDFLAGS		= -g -O -Wall -fPIC
SOFLAGS		= -shared

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= $(ROOTLIBS)  -lEG 
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard src/*.cc)
HEADERS		= $(wildcard interface/*.h)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))



all:
	$(CXX) addSystematics.cpp -o addSystematics -L${ROOTSYS}/lib $(ROOTCFLAGS) $(ROOTLIBS) -O2 -Wall 

clean:
	rm addSystematics
