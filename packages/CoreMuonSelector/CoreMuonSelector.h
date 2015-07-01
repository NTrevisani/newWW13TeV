///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                     WW ANALYSIS SELECTOR                                    /////////////
/////////////                                                                                             /////////////
/////////////                                    Nicolò Trevisani (IFCA)                                  /////////////
/////////////                                          Jun 2016                                           /////////////
/////////////                                                                                             /////////////
/////////////                              -> Adjust to a 120 width window <-                             /////////////
/////////////                                                                                             /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "PAF/computing/PAFChainItemSelector.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <vector>
#include "Riostream.h"  

const Double_t ZMASS = 91.1876;

class CoreMuonSelector: public PAFChainItemSelector{
  
 public:
  virtual ~CoreMuonSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation

  //PUWeight* fPUWeight;

  float Testing(int k);
  bool  IsTightLepton(int k);
  float MuonIsolation(int k);
  float ElectronIsolation(int k);
  bool  IsIsolatedLepton(int k);

  bool G_Debug_DefineAnalysisVariables;

  // My Declarations:
  // Define data members

  int SelectedChannel;

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
 
  std::vector<float> *std_vector_lepton_pt;      
  std::vector<float> *std_vector_jet_pt;         
  std::vector<float> *std_vector_lepton_muSIP3D; 
  std::vector<float> *std_vector_lepton_elSIP3D; 
  std::vector<float> *std_vector_lepton_id;
  std::vector<float> *std_vector_lepton_isTightMuon;
  std::vector<float> *std_vector_electron_scEta;
  std::vector<float> *std_vector_electron_deltaEtaIn;
  std::vector<float> *std_vector_electron_deltaPhiIn;
  std::vector<float> *std_vector_electron_sigmaIetaIeta;
  std::vector<float> *std_vector_electron_HoE;
  std::vector<float> *std_vector_electron_d0;
  std::vector<float> *std_vector_electron_dz;
  std::vector<float> *std_vector_electron_ooEooP;
  std::vector<float> *std_vector_electron_passConversion;
  std::vector<float> *std_vector_lepton_chargedHadronIso;
  std::vector<float> *std_vector_lepton_photonIso;
  std::vector<float> *std_vector_lepton_neutralHadronIso;
  std::vector<float> *std_vector_lepton_sumPUPt;
  std::vector<float> *std_vector_electron_effectiveArea;

  float jetRho;
  float puW;
  float effW;
  float triggW;
  float dataset;
  float baseW;
  float dphilljetjet;
  float dphilmet1;
  float dphilmet2;
  float pfType1Met;
  float pt1;
  float pt2;
  float ptll;
  float mll;
  float mth;
  float drll;
  float dphill;
  float dphilljet;
  float trkMet;
  float ch1;
  float ch2;
  float channel;
  float njet;
  float nbjet;
  float nbjettche;

  int   nextra;
  int   bveto_ip;
  int   bveto_mu;

  // VARIABLES FOR ALL EVENTS (to be initialized only once)

 
   // Histograms 
  TH1F                        *h_N_PV;
  
  // Counting histograms                                                                  
  //-------------------------------------------------------------------------

  TH1F* hWTrigger;
  TH1F* hWMetCut;
  TH1F* hWLowMinv;
  TH1F* hWZVeto;
  TH1F* hWpMetCut;
  TH1F* hWJetVeto;
  TH1F* hWnJetsBeforeBtag;
  TH1F* hWeffnJetsBeforeBtag;
  TH1F* hWnJets;
  TH1F* hWeffnJets;
  TH1F* hWnBtaggedJets;
  TH1F* hWeffnBtaggedJets;
  TH1F* hWnJetsBveto;
  TH1F* hWeffnJetsBveto;
  TH1F* hNjetsTwoLeptonsLevel;
  TH1F* hWeffnJetsBvetoAfterHt;
  TH1F* hWnJetsBvetoAfterHt;

  TH1F* hWeffTrigger;
  TH1F* hWeffMetCut;
  TH1F* hWeffLowMinv;
  TH1F* hWeffZVeto;
  TH1F* hWeffpMetCut;
  TH1F* hWeffJetVeto;
  TH1F* hWeffDeltaPhiJet;
  TH1F* hWeffSoftMuVeto;
  TH1F* hWeffExtraLepton;
  TH1F* hWeffPtll;
  TH1F* hWeffTopTagging;
 
  TH1F* hWDeltaPhiJet;
  TH1F* hWSoftMuVeto;
  TH1F* hWExtraLepton;
  TH1F* hWPtll;
  TH1F* hWTopTagging;
  
  TH1F* hPtLepton1WWLevel[4];
  TH1F* hPtLepton2WWLevel[4];
  TH1F* hPtDiLeptonWWLevel[4];
  TH1F* hMinvWWLevel[4];
  TH1F* hMtWWLevel[4];
  TH1F* hpfMetWWLevel[4];
  TH1F* hpminMetWWLevel[4];
  TH1F* hDeltaRLeptonsWWLevel[4];
  TH1F* hDeltaPhiLeptonsWWLevel[4];
  TH1F* hDPhiPtllJetWWLevel[4];
  TH1F* hSigMu[4];
  TH1F* hSigEl[4];

  TH1F* hPtLepton1WWLevelNoHt[4];
  TH1F* hPtLepton2WWLevelNoHt[4];
  TH1F* hPtDiLeptonWWLevelNoHt[4];
  TH1F* hMinvWWLevelNoHt[4];
  TH1F* hMtWWLevelNoHt[4];
  TH1F* hpfMetWWLevelNoHt[4];
  TH1F* hpminMetWWLevelNoHt[4];
  TH1F* hDeltaRLeptonsWWLevelNoHt[4];
  TH1F* hDeltaPhiLeptonsWWLevelNoHt[4];
  TH1F* hDPhiPtllJetWWLevelNoHt[4];
  TH1F* hSigMuNoHt[4];
  TH1F* hSigElNoHt[4];

  TH1F* hPtLepton1WWLevelHtPlus[4];
  TH1F* hPtLepton2WWLevelHtPlus[4];
  TH1F* hPtDiLeptonWWLevelHtPlus[4];
  TH1F* hMinvWWLevelHtPlus[4];
  TH1F* hMtWWLevelHtPlus[4];
  TH1F* hpfMetWWLevelHtPlus[4];
  TH1F* hpminMetWWLevelHtPlus[4];
  TH1F* hDeltaRLeptonsWWLevelHtPlus[4];
  TH1F* hDeltaPhiLeptonsWWLevelHtPlus[4];
  TH1F* hDPhiPtllJetWWLevelHtPlus[4];
  TH1F* hSigMuHtPlus[4];
  TH1F* hSigElHtPlus[4];

  TH1F* hHt[4];
  TH1F* hHtAfter[4];

  // TwoLeptons level histograms                                                                                                               
  //----------------------------------------------------------------------------                                                               

  TH1F* hPtLepton1TwoLeptonsLevel;
  TH1F* hPtLepton2TwoLeptonsLevel ;
  TH1F* hPtDiLeptonTwoLeptonsLevel;
  TH1F* hMinvTwoLeptonsLevel;
  TH1F* hMtTwoLeptonsLevel;
  TH1F* hpfMetTwoLeptonsLevel;
  TH1F* hpminMetTwoLeptonsLevel;
  TH1F* hDeltaRLeptonsTwoLeptonsLevel;
  TH1F* hDeltaPhiLeptonsTwoLeptonsLevel;
  TH1F* hDPhiPtllJetTwoLeptonsLevel;
  TH1F* hSigMuNoHtTwoLeptonsLevel;
  TH1F* hSigElNoHtTwoLeptonsLevel;

  TH1F *hLooseIso;

  // Input parameters
  TString                     _Signal;               // Type of Signal
  int                         _NEvents;              // Total number of events in the sample before skim
  float                       _Luminosity;           // Total luminosity
  float                       _XSection;             // Process cross section
  bool                        _IsDATA;               // True if is Data, False in case MC
  int                         _WhichRun;             // 1 in case of RunI samples. 2 In case of RunII samples.
  bool                        _Debug;                // True for verbose while debugging
  bool                        _Report;               // Count events and print final report
  TString                     _SameSign;             // Choose the type of events looking at the lepton's charge
  TString                     _FlavorChannel;        // Choose the type of events looking at the lepton's charge
  TString                     _outPath;              // Output folder

  // Weights
  float                       _factN;                // Normalization factor



 public:  
 CoreMuonSelector() : 
  
  h_N_PV(),
    
    hWTrigger(),
    hWMetCut(),
    hWLowMinv(),
    hWZVeto(),
    hWpMetCut(),
    hWJetVeto(),
    hWnJetsBeforeBtag(),
    hWeffnJetsBeforeBtag(),
    hWnJets(),
    hWeffnJets(),
    hWnBtaggedJets(),
    hWeffnBtaggedJets(),
    hWnJetsBveto(),
    hWeffnJetsBveto(),
    hNjetsTwoLeptonsLevel(),
    hWeffnJetsBvetoAfterHt(),
    
    hWDeltaPhiJet(),
    hWSoftMuVeto(),
    hWExtraLepton(),
    hWPtll(),
    hWTopTagging(),
    
    hPtLepton1WWLevel(),
    hPtLepton2WWLevel(),
    hPtDiLeptonWWLevel(),
    hMinvWWLevel(),
    hMtWWLevel(),
    hpfMetWWLevel(),
    hpminMetWWLevel(),
    hDeltaRLeptonsWWLevel(),
    hDeltaPhiLeptonsWWLevel(),
    hDPhiPtllJetWWLevel(),
    hSigMu(),
    hSigEl(),
    
    hPtLepton1WWLevelNoHt(),
    hPtLepton2WWLevelNoHt(),
    hPtDiLeptonWWLevelNoHt(),
    hMinvWWLevelNoHt(),
    hMtWWLevelNoHt(),
    hpfMetWWLevelNoHt(),
    hpminMetWWLevelNoHt(),
    hDeltaRLeptonsWWLevelNoHt(),
    hDeltaPhiLeptonsWWLevelNoHt(),
    hDPhiPtllJetWWLevelNoHt(),
    hSigMuNoHt(),
    hSigElNoHt(),
    
    hPtLepton1WWLevelHtPlus(),
    hPtLepton2WWLevelHtPlus(),
    hPtDiLeptonWWLevelHtPlus(),
    hMinvWWLevelHtPlus(),
    hMtWWLevelHtPlus(),
    hpfMetWWLevelHtPlus(),
    hpminMetWWLevelHtPlus(),
    hDeltaRLeptonsWWLevelHtPlus(),
    hDeltaPhiLeptonsWWLevelHtPlus(),
    hDPhiPtllJetWWLevelHtPlus(),
    hSigMuHtPlus(),
    hSigElHtPlus(),
    
    hHt(),
    hHtAfter(),
    
    hPtLepton1TwoLeptonsLevel(),
    hPtLepton2TwoLeptonsLevel (),
    hPtDiLeptonTwoLeptonsLevel(),
    hMinvTwoLeptonsLevel(),
    hMtTwoLeptonsLevel(),
    hpfMetTwoLeptonsLevel(),
    hpminMetTwoLeptonsLevel(),
    hDeltaRLeptonsTwoLeptonsLevel(),
    hDeltaPhiLeptonsTwoLeptonsLevel(),
    hDPhiPtllJetTwoLeptonsLevel(),
    
    hLooseIso(),
    
    _Signal(),
    _NEvents(),
    _Luminosity(),
    _XSection(),
    _IsDATA(),
    _WhichRun(),
    _Debug(),
    _Report(),
    _factN()
      { }
  
  ClassDef(CoreMuonSelector,0);
};
