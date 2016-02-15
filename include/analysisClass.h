#ifndef analysisClass_h
#define analysisClass_h

#include "baseClass.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <memory>

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

#include <boost/shared_ptr.hpp>
// typedefs
typedef boost::shared_ptr<fastjet::ClusterSequence>  ClusterSequencePtr;
typedef boost::shared_ptr<fastjet::JetDefinition>    JetDefPtr;
// For JECs
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

class analysisClass : public baseClass {
public :
  analysisClass(string * inputList, string * cutFile, string * treeName,  string *outputFileName=0, string * cutEfficFile=0);
  virtual ~analysisClass();
  void Loop();
  string fileList; //added
private :
  ClusterSequencePtr  fjClusterSeq, fjClusterSeq_shift;
  JetDefPtr           fjJetDefinition;
  // For JECs
  JetCorrectorParameters *L1Par;
  JetCorrectorParameters *L2Par;
  JetCorrectorParameters *L3Par;
  JetCorrectorParameters *L1DATAPar;
  JetCorrectorParameters *L2DATAPar;
  JetCorrectorParameters *L3DATAPar;
  JetCorrectorParameters *L2L3Residual;
  FactorizedJetCorrector *JetCorrector;
  FactorizedJetCorrector *JetCorrector_data;
  JetCorrectionUncertainty *unc;
};

#endif

#ifdef analysisClass_cxx

#endif // #ifdef analysisClass_cxx
