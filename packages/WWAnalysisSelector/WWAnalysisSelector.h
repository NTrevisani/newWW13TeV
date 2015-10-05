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

#include "/gpfs/csic_projects/cms/sw/PAF/include/PAFChainItemSelector.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <vector>
#include "Riostream.h"  

const Double_t ZMASS = 91.1876;

class WWAnalysisSelector: public PAFChainItemSelector{
  
 public:
  virtual ~WWAnalysisSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation

  //PUWeight* fPUWeight;

  float Testing(int k);
  bool  IsTightLepton(int k, TString _MuonID_);
  float MuonIsolation(int k);
  float ElectronIsolation(int k);
  bool  IsIsolatedLepton(int k);

  bool G_Debug_DefineAnalysisVariables;

  // My Declarations:
  // Define data members

  int SelectedChannel;

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
 
  std::vector<float> std_vector_lepton_pt;      
  std::vector<float> std_vector_jet_pt;         
  std::vector<float> std_vector_lepton_muSIP3D; 
  std::vector<float> std_vector_lepton_elSIP3D; 
  std::vector<float> std_vector_lepton_flavour;
  std::vector<float> std_vector_lepton_isTightMuon;
  std::vector<float> std_vector_lepton_eleIdMedium;
  std::vector<float> std_vector_lepton_chargedHadronIso;
  std::vector<float> std_vector_lepton_photonIso;
  std::vector<float> std_vector_lepton_neutralHadronIso;
  std::vector<float> std_vector_lepton_sumPUPt;
  std::vector<float> std_vector_electron_effectiveArea;
  std::vector<float> std_vector_lepton_BestTrackdxy;
  std::vector<float> std_vector_lepton_BestTrackdz;
  std::vector<float> std_vector_lepton_phi;      
  std::vector<float> std_vector_lepton_eta;      
  std::vector<float> std_vector_jet_eta;         
  std::vector<float> std_vector_jet_phi;         
  std::vector<float> std_vector_lepton_isMediumMuon;
  std::vector<float> std_vector_jet_csvv2ivf;
  std::vector<float> std_vector_jet_softMuPt;

  float GEN_weight_SM;
  float phi1;
  float phi2;
  float jetRho;
  float puW;
  float effW;
  float triggW;
  float trigger;
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
  float jetpt1;
  float jetphi1;
  float pfType1Metphi;
  float dphillmet;
  float nvtx;

  int   dphiv;
  int   nextra;
  int   bveto_ip;
  int   bveto_mu;

  // VARIABLES FOR ALL EVENTS (to be initialized only once)
  
  TTree *tree;
  
  // Histograms 
  TH1F *h_n_PV;
  
  // Counting histograms
  //----------------------------------------------------------------------------       
  
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
  
  TH1F* hWDeltaPhiJet;
  TH1F* hWSoftMuVeto;
  TH1F* hWExtraLepton;
  TH1F* hWPtll;
  TH1F* hWTopTagging;
  
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
  
  // TwoLeptons level histograms
  //----------------------------------------------------------------------------
  
  TH1F* hNVtx;
  TH1F* hPtLepton1TwoLeptonsLevel;
  TH1F* hPtLepton2TwoLeptonsLevel ;
  TH1F* hPtDiLeptonTwoLeptonsLevel;
  TH1F* hMinvTwoLeptonsLevel;
  TH1F* hMtTwoLeptonsLevel;
  TH1F* hpfMetTwoLeptonsLevel;
  TH1F* htrkMetTwoLeptonsLevel;
  TH1F* hpminMetTwoLeptonsLevel;
  TH1F* hDeltaRLeptonsTwoLeptonsLevel;
  TH1F* hDeltaPhiLeptonsTwoLeptonsLevel;
  TH1F* hDPhiPtllJetTwoLeptonsLevel;
  TH1F* hNjetsPlot1TwoLeptonsLevel;
  TH1F* hNjetsPlot2TwoLeptonsLevel;
  TH1F* hSigMuNoHtTwoLeptonsLevel;
  TH1F* hSigElNoHtTwoLeptonsLevel;
  TH1F* hDxyTwoLeptonsLevel;
  TH1F* hDzTwoLeptonsLevel;  
  
  TH1F *hLooseIso;
  TH1F *hsoftMuPt;
  TH1F *hjetPt;
  
  // WW level histograms
  //----------------------------------------------------------------------------
  
  TH1F* hWnJetsBvetoAfterHt;
  TH1F* hPtLepton1WWLevel[4];
  TH1F* hPtLepton2WWLevel[4];
  TH1F* hPtDiLeptonWWLevel[4];
  TH1F* hMinvWWLevel[4];
  TH1F* hMtWWLevel[4];
  TH1F* hpfMetWWLevel[4];
  TH1F* htrkMetWWLevel[4];
  TH1F* hpminMetWWLevel[4];
  TH1F* hDeltaRLeptonsWWLevel[4];
  TH1F* hDeltaPhiLeptonsWWLevel[4];
  TH1F* hDPhiPtllJetWWLevel[4];
  TH1F* hSigMu[4];
  TH1F* hSigEl[4];
  
  TH1F* hHt[4];
  TH1F* hHtAfter[4];
  
  //B-Veto Level Histograms
  //----------------------------------------------------------------------------    
  
  TH1F* hbVetoOld[4];
  TH1F* hbVetoMu[4];
  TH1F* hbVetoCsvv2ivfLoose[4];
  TH1F* hbVetoCsvv2ivfMedium[4];
  TH1F* hbVetoCsvv2ivfTight[4];
  TH1F* hbVetoCsvv2ivfLooseAndMu[4];
  
 public:
  
  //Additional Variables
  
  Float_t fullpmet;
  Float_t trkpmet;
  Float_t mpmet;
  Float_t Ht;
  Float_t dphijet1met;
  Float_t ratioMet;
  Float_t ptWW;
  Float_t metvar;
  Float_t njetGen;
  
  // My Declarations:OA
  // Define global variables
  
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
  TString                     _MuonID;               // ID requirement to consider a muon

  // Weights
  float                       _factN;                // Normalization factor



 public:  
 WWAnalysisSelector() : 
  
  // VARIABLES FOR ALL EVENTS (to be initialized only once)

  h_n_PV(),
  
  // Counting histograms
  //----------------------------------------------------------------------------       
  
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
  
  hWeffTrigger(),
  hWeffMetCut(),
  hWeffLowMinv(),
  hWeffZVeto(),
  hWeffpMetCut(),
  hWeffJetVeto(),
  hWeffDeltaPhiJet(),
  hWeffSoftMuVeto(),
  hWeffExtraLepton(),
  hWeffPtll(),
  hWeffTopTagging(),
  
  // TwoLeptons level histograms                                            
  //---------------------------------------------------------------------------- 
  
  hPtLepton1TwoLeptonsLevel(),
  hPtLepton2TwoLeptonsLevel (),
  hPtDiLeptonTwoLeptonsLevel(),
  hMinvTwoLeptonsLevel(),
  hMtTwoLeptonsLevel(),
  hpfMetTwoLeptonsLevel(),
  htrkMetTwoLeptonsLevel(),
  hpminMetTwoLeptonsLevel(),
  hDeltaRLeptonsTwoLeptonsLevel(),
  hDeltaPhiLeptonsTwoLeptonsLevel(),
  hDPhiPtllJetTwoLeptonsLevel(),
  hNjetsPlot1TwoLeptonsLevel(),
  hNjetsPlot2TwoLeptonsLevel(),
  hSigMuNoHtTwoLeptonsLevel(),
  hSigElNoHtTwoLeptonsLevel(),
  hDxyTwoLeptonsLevel(),
  hDzTwoLeptonsLevel(),  
  
  hLooseIso(),
  hsoftMuPt(),
  hjetPt(),

  // WW level histograms                                                         
  //---------------------------------------------------------------------------- 
  
  hPtLepton1WWLevel(),
  hPtLepton2WWLevel(),
  hPtDiLeptonWWLevel(),
  hMinvWWLevel(),
  hMtWWLevel(),
  hpfMetWWLevel(),
  htrkMetWWLevel(),
  hpminMetWWLevel(),
  hDeltaRLeptonsWWLevel(),
  hDeltaPhiLeptonsWWLevel(),
  hDPhiPtllJetWWLevel(),
  hSigMu(),
  hSigEl(),
  
  hHt(),
  hHtAfter(),
  
  //B-Veto Level Histograms
  //----------------------------------------------------------------------------    
  
  hbVetoOld(),
  hbVetoMu(),
  hbVetoCsvv2ivfLoose(),
  hbVetoCsvv2ivfMedium(),
  hbVetoCsvv2ivfTight(),
  hbVetoCsvv2ivfLooseAndMu(),
  
    _Signal(),
    _NEvents(),
    _Luminosity(),
    _XSection(),
    _IsDATA(),
    _WhichRun(),
    _Debug(),
    _Report(),
    _MuonID(),
    _factN()
      { }
  
  ClassDef(WWAnalysisSelector,0);
};
