COMP=g++ -g
FLAGS = -std=c++11 -Wno-attributes
#FLAGS += -DUSE_EXAMPLE
FLAGS += -DSAVE_ALL_HISTOGRAMS 
# FLAGS += -DCREATE_OPT_CUT_FILE
ROOTLIBS = `root-config --glibs --cflags` -lMinuit 
scram = $(shell cd ${CMSSW_BASE}; scram tool tag $(1) $(2))
INC= -I.. -I. -I./include  -I${CLHEP}/include -I${CMSSW_RELEASE_BASE}/src -I$(call scram,boost,INCLUDE) -I$(call scram,fastjet,INCLUDE)
ROOTINC= -I${ROOTSYS}/include
LIBS= -L.  ${ROOTLIBS} -L${CLHEP}/lib
SRC= ./src
SELECTIONLIB=$(SRC)/rootNtupleClass.o $(SRC)/baseClass.o $(SRC)/analysisClass.o ${CMSSW_RELEASE_BASE}/lib/${SCRAM_ARCH}/libCondFormatsJetMETObjects.so $(SRC)/jsonParser.o $(SRC)/pileupReweighter.o $(SRC)/qcdFitter.o $(SRC)/qcdFitter_V1.o  $(SRC)/likelihoodGetter.o $(SRC)/eventListHelper.o $(SRC)/utility.o $(SRC)/setTDRStyle.o #$(SRC)/CMS_lumi.o
HISTOLIB=$(SRC)/utility.o $(SRC)/setTDRStyle.o #$(SRC)/RooDCBShape.o #$(SRC)/CMS_lumi.o
ROOFITLIB= -lRooFit -lRooFitCore -lHtml -lMinuit
EXE = main
HHLIST=$(filter-out include/LinkDef.h,$(wildcard include/*.h))


# ********** TEMPLATE *************
# mainProg: mainProg.o $(SELECTIONLIB)
#	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(ROOTLIBS) -o $@  $(SELECTIONLIB) $@.o
# *********************************

all: ${EXE} histoPlotter singleHistoPlotter doPlots fit_simultaneous #CMS_lumi #superimposeDataset 
#all: histoPlotter singleHistoPlotter doPlots

main: $(SRC)/main.o $(SELECTIONLIB) 
	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(SELECTIONLIB) $(SRC)/$@.o

histoPlotter: $(SRC)/histoPlotter.o $(HISTOLIB)
	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o

singleHistoPlotter: $(SRC)/singleHistoPlotter.o $(HISTOLIB)
	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o

doPlots: $(SRC)/doPlots.o $(HISTOLIB)
	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o

#RooDCBShape: $(SRC)/RooDCBShape.o $(HISTOLIB)
#	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(ROOFITLIB) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o

fit_simultaneous: $(SRC)/fit_simultaneous.o $(HISTOLIB)
	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(ROOFITLIB) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) ./include/LinkDef.h $(SRC)/$@.o

#CMS_lumi: $(SRC)/CMS_lumi.o $(HISTOLIB)
#	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o


#superimposeDataset: $(SRC)/superimposeDataset.o $(HISTOLIB)
#	$(COMP) $(INC) $(ROOTINC) $(LIBS) $(FLAGS) `$(call scram,fastjet,FASTJET_BASE)/bin/fastjet-config --cxxflags --plugins --libs` -o $@  $(HISTOLIB) $(SRC)/$@.o

clean:
	rm -f src/*.o *.lo core core.*
	rm -f *~
	rm -f *.exe
	rm -f $(EXE)

.cpp.o:
	$(COMP) -c $(INC) $(ROOTINC) $(FLAGS) -o $@ $<

.cc.o:
	$(COMP) -c $(INC) $(ROOTINC) $(FLAGS) -o $@ $< #-m32 -c $(INC) $(ROOTINC) $(FLAGS) -o $@ $< #changed to compile new programs...

.cxx.o:
	$(COMP) -c $(INC) $(ROOTINC) $(FLAGS) -o $@ $<

.C.o:
	$(COMP) -c $(INC) $(ROOTINC) $(FLAGS) -o $@ $<


