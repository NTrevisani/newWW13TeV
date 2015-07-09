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

#include "WWAnalysisSelector.h"

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <vector>
#include "TROOT.h"
#include <iostream>
#include "TSystem.h"
#include "TLorentzVector.h"

#include "TDatabasePDG.h"

ClassImp(WWAnalysisSelector)

float dist(float phi1, float eta1, float phi2, float eta2){
  float dphi = fabs(phi1 - phi2);
  float deta = fabs(eta1 - eta2);
  if (dphi > TMath::Pi())
    dphi = 2*TMath::Pi() - dphi;
  return sqrt(dphi*dphi + deta*deta);
}


// Initialise input parameters and data members for all events
void WWAnalysisSelector::Initialise() {

  _Signal          = GetParam<TString>("Signal");
  _IsDATA          = GetParam<bool>("IsDATA");
  _NEvents         = GetParam<int>("NEvents");
  _Luminosity      = GetParam<float>("Luminosity");
  _XSection        = GetParam<float>("XSection");
  _WhichRun        = GetParam<int>("WhichRun"); 
  _Debug           = GetParam<bool>("Debug");
  _Report          = GetParam<bool>("Report");
  _SameSign        = GetParam<TString>("SameSign"); 
  _FlavorChannel   = GetParam<TString>("FlavorChannel");
  _outPath         = GetParam<TString>("OutPath");
  _MuonID          = GetParam<TString>("MuonID");

  gSystem -> Exec("mkdir " + _outPath);

  //Define weights
  _factN = 1.;
  //  if (!_IsDATA && _XSection > 0) _factN = _XSection * _Luminosity / _NEvents;

  //For counting
  if (_Report) {

  }

  SelectedChannel = -999;
  
  if      (_FlavorChannel == "MuMu") SelectedChannel =  0;
  else if (_FlavorChannel == "EE"  ) SelectedChannel =  1;
  else if (_FlavorChannel == "EMu" ) SelectedChannel =  2;
  else if (_FlavorChannel == "MuE" ) SelectedChannel =  3;
  else if (_FlavorChannel == "OF"  ) SelectedChannel =  4;
  else if (_FlavorChannel == "SF"  ) SelectedChannel =  5;
  
  else if (_FlavorChannel == "All" ) SelectedChannel = -1;
 

  //------------------------------------------------------------------------------
  // Create tree and branches
  //------------------------------------------------------------------------------
  
  tree = CreateTree("nt","nt");
  
  tree->Branch("pt1",&pt1,"pt1");
  tree->Branch("pt2",&pt2,"pt2");
  tree->Branch("pfType1Met",&pfType1Met,"pfType1Met");
  tree->Branch("trkMet",&trkMet,"trkMet");
  tree->Branch("jetpt1",&jetpt1,"jetpt1");
  tree->Branch("ptll",&ptll,"ptll");
  tree->Branch("dphill",&dphill,"dphill");
  tree->Branch("jetphi1",&jetphi1,"jetphi1");
  tree->Branch("pfType1Metphi",&pfType1Metphi,"pfType1Metphi");
  tree->Branch("mth",&mth,"mth");
  tree->Branch("mll",&mll,"mll");
  tree->Branch("mpmet",&mpmet,"mpmet");
  tree->Branch("metvar",&metvar,"metvar");
  tree->Branch("dphilljet",&dphilljet,"dphilljet");
  tree->Branch("dphillmet",&dphillmet,"dphillmet");
  tree->Branch("njet",&njet,"njet");
  tree->Branch("dphilmet1",&dphilmet1,"dphilmet1");
  tree->Branch("dphilmet2",&dphilmet2,"dphilmet2");
  tree->Branch("bveto_ip",&bveto_ip,"bveto_ip");
  tree->Branch("nbjettche",&nbjettche,"nbjettche");
  tree->Branch("dphiv",&dphiv,"dphiv");
  tree->Branch("dphijet1met",&dphijet1met,"dphijet1met");
  tree->Branch("ratioMet",&ratioMet,"ratioMet");

  tree->Branch("nvtx",&nvtx,"nvtx");
  tree->Branch("nextra",&nextra,"nextra");

  tree->Branch("std_vector_lepton_pt",&std_vector_lepton_pt,"std_vector_lepton_pt");
  tree->Branch("std_vector_jet_pt",&std_vector_jet_pt,"std_vector_jet_pt");


  //------------------------------------------------------------------------------
  // Create histos
  //------------------------------------------------------------------------------
  
  h_N_PV  = CreateH1F ("h_N_PV","h_N_PV",50,0,50); 
  h_N_PV->TH1::SetDefaultSumw2();

  // Counting histograms
  //----------------------------------------------------------------------------
  
  hWTrigger              = CreateH1F("hWTrigger",     "", 10, 0, 10);
  hWMetCut               = CreateH1F("hWMetCut",      "", 10, 0, 10);
  hWLowMinv              = CreateH1F("hWLowMinv",     "", 10, 0, 10);
  hWZVeto                = CreateH1F("hWZVeto",       "", 10, 0, 10);
  hWpMetCut              = CreateH1F("hWpMetCut",     "", 10, 0, 10);
  hWJetVeto              = CreateH1F("hWJetVeto",     "", 10, 0, 10);
  hWnJets                = CreateH1F("hWnJets",         "", 10, 0, 10);
  hWeffnJets             = CreateH1F("hWeffnJets",      "", 10, 0, 10);
  hWnBtaggedJets         = CreateH1F("hWnBtaggedJets",          "", 10, 0, 10);
  hWeffnBtaggedJets      = CreateH1F("hWeffnBtaggedJets",       "", 10, 0, 10);
  hWnJetsBveto           = CreateH1F("hWnJetsBveto",         "", 10, 0, 10);
  hWeffnJetsBveto        = CreateH1F("hWeffnJetsBveto",      "", 10, 0, 10);
  hNjetsTwoLeptonsLevel  = CreateH1F("hNjetsTwoLeptonsLevel",         "", 10, 0, 10);
  hNjetsPlot1TwoLeptonsLevel  = CreateH1F("hNjetsPlot1TwoLeptonsLevel",         "", 1000, 0, 10);
  hNjetsPlot2TwoLeptonsLevel  = CreateH1F("hNjetsPlot2TwoLeptonsLevel",         "", 1000, 0, 10);
  hWeffnJetsBvetoAfterHt = CreateH1F("hWeffnJetsBvetoAfterHt",      "", 10, 0, 10);

  hWDeltaPhiJet          = CreateH1F("hWDeltaPhiJet", "", 10, 0, 10);
  hWSoftMuVeto           = CreateH1F("hWSoftMuVeto",  "", 10, 0, 10);
  hWExtraLepton          = CreateH1F("hWExtraLepton", "", 10, 0, 10);
  hWPtll                 = CreateH1F("hWPtll",        "", 10, 0, 10);
  hWTopTagging           = CreateH1F("hWTopTagging",  "", 10, 0, 10);

  hWeffTrigger           = CreateH1F("hWeffTrigger",     "", 10, 0, 10);
  hWeffMetCut            = CreateH1F("hWeffMetCut",      "", 10, 0, 10);
  hWeffLowMinv           = CreateH1F("hWeffLowMinv",     "", 10, 0, 10);
  hWeffZVeto             = CreateH1F("hWeffZVeto",       "", 10, 0, 10);
  hWeffpMetCut           = CreateH1F("hWeffpMetCut",     "", 10, 0, 10);
  hWeffJetVeto           = CreateH1F("hWeffJetVeto",     "", 10, 0, 10);
  hWeffDeltaPhiJet       = CreateH1F("hWeffDeltaPhiJet", "", 10, 0, 10);
  hWeffSoftMuVeto        = CreateH1F("hWeffSoftMuVeto",  "", 10, 0, 10);
  hWeffExtraLepton       = CreateH1F("hWeffExtraLepton", "", 10, 0, 10);
  hWeffPtll              = CreateH1F("hWeffPtll",        "", 10, 0, 10);
  hWeffTopTagging        = CreateH1F("hWeffTopTagging",  "", 10, 0, 10);

  hLooseIso              = CreateH1F("hLooseIso",        "", 100, 0, 10);

  // WW level histograms
  //---------------------------------------------------------------------------

  for (Int_t nC=0; nC<4; nC++) {

    hPtLepton1WWLevel[nC]             = CreateH1F(Form("hPtLepton1WWLevel%.1i", nC),             "", 1000, 0, 1000);
    hPtLepton2WWLevel[nC]             = CreateH1F(Form("hPtLepton2WWLevel%.1i", nC),             "", 1000, 0, 1000);
    hPtDiLeptonWWLevel[nC]            = CreateH1F(Form("hPtDiLeptonWWLevel%.1i", nC),            "", 1000, 0, 1000);
    hMinvWWLevel[nC]                  = CreateH1F(Form("hMinvWWLevel%.1i", nC),                  "", 1000, 0, 1000);
    hMtWWLevel[nC]                    = CreateH1F(Form("hMtWWLevel%.1i", nC),                    "", 1000, 0, 1000);
    hpfMetWWLevel[nC]                 = CreateH1F(Form("hpfMetWWLevel%.1i", nC),                 "", 1000, 0, 1000);
    hpminMetWWLevel[nC]               = CreateH1F(Form("hpminMetWWLevel%.1i", nC),               "", 1000, 0, 1000);
    hDeltaRLeptonsWWLevel[nC]         = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),         "",   50, 0,    5);
    hDeltaPhiLeptonsWWLevel[nC]       = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC),       "",   32, 0,  3.2);
    hDPhiPtllJetWWLevel[nC]           = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),           "",   32, 0,  3.2);
    hSigEl[nC]                        = CreateH1F(Form("hSigEl%.1i", nC),                        "", 1000, 0, 1000);
    hSigMu[nC]                        = CreateH1F(Form("hSigMu%.1i", nC),                        "", 1000, 0, 1000);

    hPtLepton1WWLevelNoHt[nC]         = CreateH1F(Form("hPtLepton1WWLevelNoHt%.1i", nC),         "", 1000, 0, 1000);
    hPtLepton2WWLevelNoHt[nC]         = CreateH1F(Form("hPtLepton2WWLevelNoHt%.1i", nC),         "", 1000, 0, 1000);
    hPtDiLeptonWWLevelNoHt[nC]        = CreateH1F(Form("hPtDiLeptonWWLevelNoHt%.1i", nC),        "", 1000, 0, 1000);
    hMinvWWLevelNoHt[nC]              = CreateH1F(Form("hMinvWWLevelNoHt%.1i", nC),              "", 1000, 0, 1000);
    hMtWWLevelNoHt[nC]                = CreateH1F(Form("hMtWWLevelNoHt%.1i", nC),                "", 1000, 0, 1000);
    hpfMetWWLevelNoHt[nC]             = CreateH1F(Form("hpfMetWWLevelNoHt%.1i", nC),             "", 1000, 0, 1000);
    hpminMetWWLevelNoHt[nC]           = CreateH1F(Form("hpminMetWWLevelNoHt%.1i", nC),           "", 1000, 0, 1000);
    hDeltaRLeptonsWWLevelNoHt[nC]     = CreateH1F(Form("hDeltaRLeptonsWWLevelNoHt%.1i", nC),     "",   50, 0,    5);
    hDeltaPhiLeptonsWWLevelNoHt[nC]   = CreateH1F(Form("hDeltaPhiLeptonsWWLevelNoHt%.1i", nC),   "",   32, 0,  3.2);
    hDPhiPtllJetWWLevelNoHt[nC]       = CreateH1F(Form("hDPhiPtllJetWWLevelNoHt%.1i", nC),       "",   32, 0,  3.2);
    hSigElNoHt[nC]                    = CreateH1F(Form("hSigElNoHt%.1i", nC),                    "", 1000, 0, 1000);
    hSigMuNoHt[nC]                    = CreateH1F(Form("hSigMuNoHt%.1i", nC),                    "", 1000, 0, 1000);

    hPtLepton1WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton1WWLevelHtPlus%.1i", nC),       "", 1000, 0, 1000);
    hPtLepton2WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton2WWLevelHtPlus%.1i", nC),       "", 1000, 0, 1000);
    hPtDiLeptonWWLevelHtPlus[nC]      = CreateH1F(Form("hPtDiLeptonWWLevelHtPlus%.1i", nC),      "", 1000, 0, 1000);
    hMinvWWLevelHtPlus[nC]            = CreateH1F(Form("hMinvWWLevelHtPlus%.1i", nC),            "", 1000, 0, 1000);
    hMtWWLevelHtPlus[nC]              = CreateH1F(Form("hMtWWLevelHtPlus%.1i", nC),              "", 1000, 0, 1000);
    hpfMetWWLevelHtPlus[nC]           = CreateH1F(Form("hpfMetWWLevelHtPlus%.1i", nC),           "", 1000, 0, 1000);
    hpminMetWWLevelHtPlus[nC]         = CreateH1F(Form("hpminMetWWLevelHtPlus%.1i", nC),         "", 1000, 0, 1000);
    hDeltaRLeptonsWWLevelHtPlus[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevelHtPlus%.1i", nC),   "",   50, 0,    5);
    hDeltaPhiLeptonsWWLevelHtPlus[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevelHtPlus%.1i", nC), "",   32, 0,  3.2);
    hDPhiPtllJetWWLevelHtPlus[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevelHtPlus%.1i", nC),     "",   32, 0,  3.2);
    hSigElHtPlus[nC]                  = CreateH1F(Form("hSigElHtPlus%.1i", nC),                  "", 1000, 0, 1000);
    hSigMuHtPlus[nC]                  = CreateH1F(Form("hSigMuHtPlus%.1i", nC),                  "", 1000, 0, 1000);

    hHt[nC]                           = CreateH1F(Form("hHt%.1i",               nC),             "", 3000, 0, 3000);
    hHtAfter[nC]                      = CreateH1F(Form("hHtAfter%.1i",          nC),             "", 3000, 0, 3000);
  }

  hWnJetsBvetoAfterHt                 = CreateH1F("hWnJetsBvetoAfterHt","",10,0,10);

  // TwoLeptons level histograms
  //----------------------------------------------------------------------------

  hPtLepton1TwoLeptonsLevel           = CreateH1F("hPtLepton1TwoLeptonsLevel",                   "", 1000, 0, 1000);
  hPtLepton2TwoLeptonsLevel           = CreateH1F("hPtLepton2TwoLeptonsLevel",                   "", 1000, 0, 1000);
  hPtDiLeptonTwoLeptonsLevel          = CreateH1F("hPtDiLeptonTwoLeptonsLevel",                  "", 1000, 0, 1000);
  hMinvTwoLeptonsLevel                = CreateH1F("hMinvTwoLeptonsLevel",                        "", 1000, 0, 1000);
  hMtTwoLeptonsLevel                  = CreateH1F("hMtTwoLeptonsLevel",                          "", 1000, 0, 1000);
  hpfMetTwoLeptonsLevel               = CreateH1F("hpfMetTwoLeptonsLevel",                       "", 1000, 0, 1000);
  hpminMetTwoLeptonsLevel             = CreateH1F("hpminMetTwoLeptonsLevel",                     "", 1000, 0, 1000);
  hDeltaRLeptonsTwoLeptonsLevel       = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",               "",   50, 0,    5);
  hDeltaPhiLeptonsTwoLeptonsLevel     = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel",             "",   32, 0,  3.2);
  hDPhiPtllJetTwoLeptonsLevel         = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",                 "",   32, 0,  3.2);
  hSigMuNoHtTwoLeptonsLevel           = CreateH1F("hSigMuNoHtTwoLeptonsLevel",                   "", 1000, 0, 1000);
  hSigElNoHtTwoLeptonsLevel           = CreateH1F("hSigElNoHtTwoLeptonsLevel",                   "", 1000, 0, 1000);

  //Vectors
  //----------------------------------------------------------------------------

  std_vector_lepton_pt               = new std::vector<float> ();      
  std_vector_jet_pt                  = new std::vector<float> ();       
  std_vector_lepton_muSIP3D          = new std::vector<float> (); 
  std_vector_lepton_elSIP3D          = new std::vector<float> ();
  std_vector_lepton_id               = new std::vector<float> ();
  std_vector_lepton_isTightMuon      = new std::vector<float> ();
  std_vector_electron_scEta          = new std::vector<float> ();
  std_vector_electron_deltaEtaIn     = new std::vector<float> ();
  std_vector_electron_deltaPhiIn     = new std::vector<float> ();
  std_vector_electron_sigmaIetaIeta  = new std::vector<float> ();
  std_vector_electron_HoE            = new std::vector<float> ();
  std_vector_electron_d0             = new std::vector<float> ();
  std_vector_electron_dz             = new std::vector<float> ();
  std_vector_electron_ooEooP         = new std::vector<float> ();
  std_vector_electron_passConversion = new std::vector<float> ();
  std_vector_electron_deltaEtaIn     = new std::vector<float> ();
  std_vector_electron_deltaPhiIn     = new std::vector<float> ();
  std_vector_electron_sigmaIetaIeta  = new std::vector<float> ();
  std_vector_electron_HoE            = new std::vector<float> ();
  std_vector_electron_d0             = new std::vector<float> ();
  std_vector_electron_dz             = new std::vector<float> ();
  std_vector_electron_ooEooP         = new std::vector<float> ();
  std_vector_electron_passConversion = new std::vector<float> ();
  std_vector_lepton_chargedHadronIso = new std::vector<float> ();
  std_vector_lepton_photonIso        = new std::vector<float> ();
  std_vector_lepton_neutralHadronIso = new std::vector<float> ();
  std_vector_lepton_sumPUPt          = new std::vector<float> ();
  std_vector_lepton_chargedHadronIso = new std::vector<float> ();
  std_vector_lepton_photonIso        = new std::vector<float> ();
  std_vector_lepton_neutralHadronIso = new std::vector<float> ();
  std_vector_electron_effectiveArea  = new std::vector<float> ();
  std_vector_lepton_BestTrackdxy     = new std::vector<float> ();
  std_vector_lepton_BestTrackdz      = new std::vector<float> ();
  std_vector_lepton_phi              = new std::vector<float> ();
  std_vector_lepton_eta              = new std::vector<float> ();
  std_vector_jet_eta                 = new std::vector<float> ();
  std_vector_jet_phi                 = new std::vector<float> ();
  std_vector_lepton_isMediumMuon     = new std::vector<float> ();
}

void WWAnalysisSelector::InsideLoop() {

  //Assigning values to vectors
  //----------------------------------------------------------------------------

  for(int i = 0; i < GetSizeOf("std_vector_lepton_pt"); ++i)
    std_vector_lepton_pt -> push_back( Get<float>("std_vector_lepton_pt", i) );

  for(int i = 0; i < GetSizeOf("std_vector_jet_pt"); ++i)
    std_vector_jet_pt -> push_back( Get<float>("std_vector_jet_pt",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_muSIP3D"); ++i)
    std_vector_lepton_muSIP3D -> push_back( Get<float>("std_vector_lepton_muSIP3D",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_elSIP3D"); ++i)
    std_vector_lepton_elSIP3D -> push_back( Get<float>("std_vector_lepton_elSIP3D",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_HoE"); ++i)
    std_vector_electron_HoE -> push_back( Get<float>("std_vector_electron_HoE",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_deltaEtaIn"); ++i)
    std_vector_electron_deltaEtaIn -> push_back( Get<float>("std_vector_electron_deltaEtaIn",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_d0"); ++i)
    std_vector_electron_d0 -> push_back( Get<float>("std_vector_electron_d0",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_sigmaIetaIeta"); ++i)
    std_vector_electron_sigmaIetaIeta -> push_back( Get<float>("std_vector_electron_sigmaIetaIeta",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_dz"); ++i)
    std_vector_electron_dz -> push_back( Get<float>("std_vector_electron_dz",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_ooEooP"); ++i)
    std_vector_electron_ooEooP -> push_back( Get<float>("std_vector_electron_ooEooP",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_ooEooP"); ++i)
    std_vector_electron_ooEooP -> push_back( Get<float>("std_vector_electron_ooEooP",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_deltaPhiIn"); ++i)
    std_vector_electron_deltaPhiIn -> push_back( Get<float>("std_vector_electron_deltaPhiIn",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_passConversion"); ++i)
    std_vector_electron_passConversion -> push_back( Get<float>("std_vector_electron_passConversion",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_dz"); ++i)
    std_vector_electron_dz -> push_back( Get<float>("std_vector_electron_dz",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_d0"); ++i)
    std_vector_electron_d0 -> push_back( Get<float>("std_vector_electron_d0",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_HoE"); ++i)
    std_vector_electron_HoE -> push_back( Get<float>("std_vector_electron_HoE",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_sigmaIetaIeta"); ++i)
    std_vector_electron_sigmaIetaIeta -> push_back( Get<float>("std_vector_electron_sigmaIetaIeta",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_deltaPhiIn"); ++i)
    std_vector_electron_deltaPhiIn -> push_back( Get<float>("std_vector_electron_deltaPhiIn",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_deltaEtaIn"); ++i)
    std_vector_electron_deltaEtaIn -> push_back( Get<float>("std_vector_electron_deltaEtaIn",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_scEta"); ++i)
    std_vector_electron_scEta -> push_back( Get<float>("std_vector_electron_scEta",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_isTightMuon"); ++i)
    std_vector_lepton_isTightMuon -> push_back( Get<float>("std_vector_lepton_isTightMuon",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_id"); ++i)
    std_vector_lepton_id -> push_back( Get<float>("std_vector_lepton_id",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_passConversion"); ++i)
    std_vector_electron_passConversion -> push_back( Get<float>("std_vector_electron_passConversion",i) );

  for(int i = 0; i < GetSizeOf("std_vector_electron_effectiveArea"); ++i)
    std_vector_electron_effectiveArea -> push_back( Get<float>("std_vector_electron_effectiveArea",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_neutralHadronIso"); ++i)
    std_vector_lepton_neutralHadronIso -> push_back( Get<float>("std_vector_lepton_neutralHadronIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_photonIso"); ++i)
    std_vector_lepton_photonIso -> push_back( Get<float>("std_vector_lepton_photonIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_chargedHadronIso"); ++i)
    std_vector_lepton_chargedHadronIso -> push_back( Get<float>("std_vector_lepton_chargedHadronIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_sumPUPt"); ++i)
    std_vector_lepton_sumPUPt -> push_back( Get<float>("std_vector_lepton_sumPUPt",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_neutralHadronIso"); ++i)
    std_vector_lepton_neutralHadronIso -> push_back( Get<float>("std_vector_lepton_neutralHadronIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_photonIso"); ++i)
    std_vector_lepton_photonIso -> push_back( Get<float>("std_vector_lepton_photonIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_chargedHadronIso"); ++i)
    std_vector_lepton_chargedHadronIso -> push_back( Get<float>("std_vector_lepton_chargedHadronIso",i) );

  for(int i = 0; i < GetSizeOf("std_vector_lepton_BestTrackdxy"); ++i)
    std_vector_lepton_BestTrackdxy -> push_back( Get<float>("std_vector_lepton_BestTrackdxy",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_lepton_BestTrackdz"); ++i)
    std_vector_lepton_BestTrackdz -> push_back( Get<float>("std_vector_lepton_BestTrackdz",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_lepton_phi"); ++i)
    std_vector_lepton_phi -> push_back( Get<float>("std_vector_lepton_phi",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_lepton_eta"); ++i)
    std_vector_lepton_eta -> push_back( Get<float>("std_vector_lepton_eta",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_jet_eta"); ++i)
    std_vector_jet_eta -> push_back( Get<float>("std_vector_jet_eta",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_jet_phi"); ++i)
    std_vector_jet_phi -> push_back( Get<float>("std_vector_jet_phi",i) );     

  for(int i = 0; i < GetSizeOf("std_vector_lepton_isMediumMuon"); ++i)
    std_vector_lepton_isMediumMuon -> push_back( Get<float>("std_vector_lepton_isMediumMuon",i) );

  //Assigning values to float variables
  //----------------------------------------------------------------------------

  Assign("puW",puW);  
  Assign("effW",effW);
  Assign("triggW",triggW);
  Assign("dataset",dataset);
  Assign("baseW",baseW);
  Assign("dphilljetjet",dphilljetjet);
  Assign("dphilmet1",dphilmet1);
  Assign("dphilmet2",dphilmet2);
  Assign("pfType1Met",pfType1Met);  
  Assign("pt1",pt1);
  Assign("pt2",pt2);
  Assign("ptll",ptll);
  Assign("mll",mll);
  Assign("mth",mth);
  Assign("drll",drll);
  Assign("dphill",dphill);
  Assign("dphilljet",dphilljet);
  Assign("trkMet",trkMet);
  Assign("ch1",ch1);
  Assign("ch2",ch2);
  Assign("jetRho",jetRho);
  Assign("channel",channel);
  Assign("njet",njet);
  Assign("nbjet",nbjet);
  Assign("nbjettche",nbjettche);
  Assign("jetpt1",jetpt1);
  Assign("jetphi1",jetphi1);
  Assign("pfType1Metphi",pfType1Metphi);
  Assign("dphillmet",dphillmet);

  //Assigning values to integer variables
  //----------------------------------------------------------------------------
  Assign("nvtx",nvtx);
  Assign("bveto_ip",bveto_ip);
  Assign("bveto_mu",bveto_mu);

  //Creating the variables we need
  //----------------------------------------------------------------------------

  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * _Luminosity;
  
  h_N_PV -> Fill(1,efficiencyW);  
  
  dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  int jetbin = njet;
  
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  Float_t fullpmet = 0;
  Float_t trkpmet = 0;
  
  if (dphimin < TMath::Pi() / 2)
    fullpmet = pfType1Met * sin(dphimin);
  else 
    fullpmet = pfType1Met;

  if (dphimin < TMath::Pi() / 2)
    trkpmet = trkMet * sin(dphimin);
  else
    trkpmet = trkMet;
  
   mpmet = min(trkpmet,fullpmet);

   metvar = (njet <= 1) ? mpmet : pfType1Met;
  
   //building Ht
   Float_t Ht = 0.;
   Ht = std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;
   
   if(njet > 10) njet = 10;
   for (int i = 0; i < njet; ++i)
     if(std_vector_jet_pt->at(i) > 0)
       Ht += std_vector_jet_pt->at(i);

   //defining nextra
   nextra = 0;
   for (int i = 0; i < std_vector_lepton_pt -> size(); ++i)
     if(std_vector_lepton_pt -> at(i) > 10)
       ++nextra;
   nextra = nextra -2;
   
   //building dRjet1
   Float_t dRjet1  = 100.;
   TLorentzVector vjet1(0.,0.,0.,0.);
   TLorentzVector vlep1(0.,0.,0.,0.);
   if(std_vector_lepton_pt->at(0) > 0.){
     vlep1.SetPtEtaPhiM(std_vector_lepton_pt->at(0),std_vector_lepton_eta->at(0),std_vector_lepton_phi->at(0),0.);
     for (int i = 0; i < std_vector_jet_pt -> size(); ++i)
       if(std_vector_jet_pt -> at(i) > 0.){
	 vjet1.SetPtEtaPhiM(std_vector_jet_pt->at(i),std_vector_jet_eta->at(i),std_vector_jet_phi->at(i),0.);
	 //cout<<i<<":"<<vlep1.DeltaR(vjet1)<<","<<dist(std_vector_lepton_phi->at(0),std_vector_lepton_eta->at(0),std_vector_jet_phi->at(i),std_vector_jet_eta->at(i))<<endl;
	 if( vlep1.DeltaR(vjet1) < dRjet1){
	   dRjet1 = vlep1.DeltaR(vjet1);
	 }
       }
   }

   //building dRjet2
   Float_t dRjet2  = 100.;
   TLorentzVector vjet2(0.,0.,0.,0.);
   TLorentzVector vlep2(0.,0.,0.,0.);
   if(std_vector_lepton_pt->at(1) > 0.){
     vlep2.SetPtEtaPhiM(std_vector_lepton_pt->at(1),std_vector_lepton_phi->at(1),std_vector_lepton_eta->at(1),0.);
     for (int i = 0; i < njet; ++i)
       if(std_vector_jet_pt -> at(i) > 0.){
	 vjet2.SetPtEtaPhiM(std_vector_jet_pt->at(i),std_vector_jet_eta->at(i),std_vector_jet_phi->at(i),0.);
	 if( vlep2.DeltaR(vjet2) < dRjet2){
	   dRjet2 = vlep2.DeltaR(vjet2);
	 }
       }
   }
   
   //Building dphijet1met
   dphijet1met = 0.;
   if (jetphi1 > 0 && pfType1Metphi > 0){
     dphijet1met = fabs(jetphi1 - pfType1Metphi);
     if (dphijet1met > TMath::Pi()) dphijet1met = 2*TMath::Pi() - dphijet1met;
   }

   //Building RatioMet   
   ratioMet = 0.;
   if (pfType1Met > 0 && trkMet > 0)
     ratioMet = pfType1Met / sqrt(pfType1Met + trkMet);


   // The selection begins here
   //--------------------------------------------------------------------------
   if (std_vector_lepton_pt->at(0) > 20)
     if (std_vector_lepton_pt->at(1) > 20) 
       if ((_SameSign == "SS" && ch1*ch2 > 0) || (_SameSign == "OS" && ch1*ch2 < 0))
	 if ( (SelectedChannel == -1)                                     || 
	      (channel == SelectedChannel)                                || 
	      (SelectedChannel == 4 && (channel == 2 || channel == 3) )   || 
	      (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
	      ){
  
	   if (IsTightLepton(0) && !IsTightLepton(1))
	     hLooseIso -> Fill(ElectronIsolation(1), totalW);
	   if (IsTightLepton(1) && !IsTightLepton(0))
	     hLooseIso -> Fill(ElectronIsolation(0), totalW);
	   
	   if (IsIsolatedLepton(0))
	     if (IsIsolatedLepton(1))
	       if (IsTightLepton(0))
		 if (IsTightLepton(1)){
	       
		   tree->Fill();
	       //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	       //
	       // Main analisis
	       //
	       //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	       
	       hWTrigger   ->Fill(1, totalW); 
	       hWeffTrigger->Fill(1, efficiencyW);
	       
	       hPtLepton1TwoLeptonsLevel      ->Fill(pt1,        totalW);
	       hPtLepton2TwoLeptonsLevel      ->Fill(pt2,        totalW);
	       hPtDiLeptonTwoLeptonsLevel     ->Fill(ptll,       totalW);
	       hMinvTwoLeptonsLevel           ->Fill(mll,        totalW);
	       hMtTwoLeptonsLevel             ->Fill(mth,        totalW);
	       hpfMetTwoLeptonsLevel          ->Fill(pfType1Met, totalW);
	       hpminMetTwoLeptonsLevel        ->Fill(mpmet,      totalW);
	       hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,       totalW);
	       hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,     totalW);
	       hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet,  totalW);
	       hNjetsTwoLeptonsLevel          ->Fill(njet,       totalW);
	       hNjetsPlot1TwoLeptonsLevel     ->Fill(dRjet1,     totalW);
	       hNjetsPlot2TwoLeptonsLevel     ->Fill(dRjet2,     totalW);
	       hSigMuNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
	       hSigElNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);

	       if (nextra == 0) {
		 
		 hWExtraLepton->Fill(1, totalW);
		 hWeffExtraLepton->Fill(1, efficiencyW);

		 if (pfType1Met > 20 ) { // removed for differential xsec
		   
		   hWMetCut->Fill(1, totalW);
		   hWeffMetCut->Fill(1, efficiencyW);
			       
		   if (mll > 12) {
		     
		     hWLowMinv->Fill(1, totalW);
		     hWeffLowMinv->Fill(1, efficiencyW);
			       
		     //zveto (in case of same flavour)
		     if ( (fabs(ZMASS - mll) > 15 && 
			   ( metvar > 45 ) )      ||
			  channel == 2            ||
			  channel == 3            ){
			    
		       hWZVeto->Fill(1, totalW);
		       hWeffZVeto->Fill(1, efficiencyW);
		       
		       if (mpmet > 20){
			 
			 hWpMetCut->Fill(1, totalW);
			 hWeffpMetCut->Fill(1, efficiencyW);

			 if (dphiv || channel == 2 || channel == 3) {
			   
			   hWDeltaPhiJet->Fill(1, totalW);
			   hWeffDeltaPhiJet->Fill(1, efficiencyW);

			   if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) {
			     
			     hWPtll->Fill(1, totalW);			    
			     hWeffPtll->Fill(1, efficiencyW);			    
			     
			     hWnJets->Fill(njet, totalW);
			     hWeffnJets->Fill(njet, efficiencyW);
			     
			     hWnBtaggedJets->Fill(nbjet, totalW);
			     hWeffnBtaggedJets->Fill(nbjet, efficiencyW);
			     
			     hHt[3]->Fill(Ht,totalW);				    

			     for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			       if (jetbin >= 3) jetbin = 2;
			       if(jetNumber == jetbin){
				 hHt[jetNumber]->Fill(Ht,totalW);				    
			       }
			     }
			     
			     //b-veto
			     if (bveto_ip == 1 && nbjettche == 0) {
			       
			       hWTopTagging->Fill(1, totalW);
			       hWeffTopTagging->Fill(1, efficiencyW);
			       hHt[2]->Fill(Ht,totalW);			       

			       //b-veto (now not operative)
			       /*if (bveto_mu == 1) */{
				 
				 hWSoftMuVeto->Fill(1, totalW);
				 hWeffSoftMuVeto->Fill(1, efficiencyW);
				 
				 hHtAfter[3]->Fill(Ht,totalW);				    
				 
				 hPtLepton1WWLevelNoHt[3]      ->Fill(pt1,       totalW);
				 hPtLepton2WWLevelNoHt[3]      ->Fill(pt2,       totalW);
				 hPtDiLeptonWWLevelNoHt[3]     ->Fill(ptll,      totalW);
				 hMinvWWLevelNoHt[3]           ->Fill(mll,       totalW);
				 hMtWWLevelNoHt[3]             ->Fill(mth,       totalW);
				 hpfMetWWLevelNoHt[3]          ->Fill(pfType1Met,totalW);
				 hpminMetWWLevelNoHt[3]        ->Fill(mpmet,     totalW);
				 hDeltaRLeptonsWWLevelNoHt[3]  ->Fill(drll,      totalW);
				 hDeltaPhiLeptonsWWLevelNoHt[3]->Fill(dphill,    totalW);
				 hDPhiPtllJetWWLevelNoHt[3]    ->Fill(dphilljet, totalW);
				 hWnJetsBveto                  ->Fill(njet,      totalW);
				 hWnJetsBveto                  ->Fill(njet, efficiencyW);
				 hSigMuNoHt[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				 hSigElNoHt[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);

				 //bveto Ht 
				 if(Ht < 218){
				   hPtLepton1WWLevel[3]      ->Fill(pt1,       totalW);
				   hPtLepton2WWLevel[3]      ->Fill(pt2,       totalW);
				   hPtDiLeptonWWLevel[3]     ->Fill(ptll,      totalW);
				   hMinvWWLevel[3]           ->Fill(mll,       totalW);
				   hMtWWLevel[3]             ->Fill(mth,       totalW);
				   hpfMetWWLevel[3]          ->Fill(pfType1Met,totalW);
				   hpminMetWWLevel[3]        ->Fill(mpmet,     totalW);
				   hDeltaRLeptonsWWLevel[3]  ->Fill(drll,      totalW);
				   hDeltaPhiLeptonsWWLevel[3]->Fill(dphill,    totalW);
				   hDPhiPtllJetWWLevel[3]    ->Fill(dphilljet, totalW);
				   hWnJetsBvetoAfterHt       ->Fill(njet, efficiencyW);					
				   hSigMu[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				   hSigEl[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				 }
				 
				 //bveto Ht 
				 if(Ht > 218){
				   
				   hPtLepton1WWLevelHtPlus[3]      ->Fill(pt1,       totalW);
				   hPtLepton2WWLevelHtPlus[3]      ->Fill(pt2,       totalW);
				   hPtDiLeptonWWLevelHtPlus[3]     ->Fill(ptll,      totalW);
				   hMinvWWLevelHtPlus[3]           ->Fill(mll,       totalW);
				   hMtWWLevelHtPlus[3]             ->Fill(mth,       totalW);
				   hpfMetWWLevelHtPlus[3]          ->Fill(pfType1Met,totalW);
				   hpminMetWWLevelHtPlus[3]        ->Fill(mpmet,     totalW);
				   hDeltaRLeptonsWWLevelHtPlus[3]  ->Fill(drll,      totalW);
				   hDeltaPhiLeptonsWWLevelHtPlus[3]->Fill(dphill,    totalW);
				   hDPhiPtllJetWWLevelHtPlus[3]    ->Fill(dphilljet, totalW);
				   hSigMuHtPlus[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				   hSigElHtPlus[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				 }
				 
				 for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
				   if (jetbin >= 3) jetbin = 2;
				   if(jetNumber == jetbin){
				     
				     hHtAfter[jetNumber]                   ->Fill(Ht,        totalW);				    
				     hPtLepton1WWLevelNoHt[jetNumber]      ->Fill(pt1,       totalW);
				     hPtLepton2WWLevelNoHt[jetNumber]      ->Fill(pt2,       totalW);
				     hPtDiLeptonWWLevelNoHt[jetNumber]     ->Fill(ptll,      totalW);
				     hMinvWWLevelNoHt[jetNumber]           ->Fill(mll,       totalW);
				     hMtWWLevelNoHt[jetNumber]             ->Fill(mth,       totalW);
				     hpfMetWWLevelNoHt[jetNumber]          ->Fill(pfType1Met,totalW);
				     hpminMetWWLevelNoHt[jetNumber]        ->Fill(mpmet,     totalW);
				     hDeltaRLeptonsWWLevelNoHt[jetNumber]  ->Fill(drll,      totalW);
				     hDeltaPhiLeptonsWWLevelNoHt[jetNumber]->Fill(dphill,    totalW);
				     hDPhiPtllJetWWLevelNoHt[jetNumber]    ->Fill(dphilljet, totalW);
				     hSigMuNoHt[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				     hSigElNoHt[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				     
				     //bveto Ht  
				     if(Ht < 218){
				       
				       hPtLepton1WWLevel[jetNumber]      ->Fill(pt1,       totalW);
				       hPtLepton2WWLevel[jetNumber]      ->Fill(pt2,       totalW);
				       hPtDiLeptonWWLevel[jetNumber]     ->Fill(ptll,      totalW);
				       hMinvWWLevel[jetNumber]           ->Fill(mll,       totalW);
				       hMtWWLevel[jetNumber]             ->Fill(mth,       totalW);
				       hpfMetWWLevel[jetNumber]          ->Fill(pfType1Met,totalW);
				       hpminMetWWLevel[jetNumber]        ->Fill(mpmet,     totalW);
				       hDeltaRLeptonsWWLevel[jetNumber]  ->Fill(drll,      totalW);
				       hDeltaPhiLeptonsWWLevel[jetNumber]->Fill(dphill,    totalW);
				       hDPhiPtllJetWWLevel[jetNumber]    ->Fill(dphilljet, totalW);
				       hSigMu[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				       hSigEl[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				     }
				     
				     //bveto Ht  
				     if(Ht > 218){
				       
				       hPtLepton1WWLevelHtPlus[jetNumber]      ->Fill(pt1,       totalW);
				       hPtLepton2WWLevelHtPlus[jetNumber]      ->Fill(pt2,       totalW);
				       hPtDiLeptonWWLevelHtPlus[jetNumber]     ->Fill(ptll,      totalW);
				       hMinvWWLevelHtPlus[jetNumber]           ->Fill(mll,       totalW);
				       hMtWWLevelHtPlus[jetNumber]             ->Fill(mth,       totalW);
				       hpfMetWWLevelHtPlus[jetNumber]          ->Fill(pfType1Met,totalW);
				       hpminMetWWLevelHtPlus[jetNumber]        ->Fill(mpmet,     totalW);
				       hDeltaRLeptonsWWLevelHtPlus[jetNumber]  ->Fill(drll,      totalW);
				       hDeltaPhiLeptonsWWLevelHtPlus[jetNumber]->Fill(dphill,    totalW);
				       hDPhiPtllJetWWLevelHtPlus[jetNumber]    ->Fill(dphilljet, totalW);
				       hSigMuHtPlus[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				       hSigElHtPlus[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				     }
				   }  					
				 }
			       }
			     }
			   }
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	 }
   
   // Define Normalization Factor for MC samples 
   //------------------------------------------------------------------------------
   
   // Define weights
   //------------------------------------------------------------------------------
   
   float pileupweight = 1;
   
   //  if (!IsDATA)
   // pileupweight = fPUWeight->GetWeight(T_Event_nPU);
   
   double factN = 1;
   
   if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents;
   
   //factN = factN*pileupweight;
   //------------------------------------------------------------------------------
   
   // Init variables
   //------------------------------------------------------------------------------

   std_vector_lepton_pt               -> clear();
   std_vector_jet_pt                  -> clear();       
   std_vector_lepton_muSIP3D          -> clear(); 
   std_vector_lepton_elSIP3D          -> clear();
   std_vector_lepton_id               -> clear();
   std_vector_lepton_isTightMuon      -> clear();
   std_vector_electron_scEta          -> clear();          
   std_vector_electron_deltaEtaIn     -> clear();     
   std_vector_electron_deltaPhiIn     -> clear();     
   std_vector_electron_sigmaIetaIeta  -> clear();  
   std_vector_electron_HoE            -> clear();            
   std_vector_electron_d0             -> clear();             
   std_vector_electron_dz             -> clear();             
   std_vector_electron_ooEooP         -> clear();         
   std_vector_electron_passConversion -> clear(); 
   std_vector_electron_deltaEtaIn     -> clear();     
   std_vector_electron_deltaPhiIn     -> clear();     
   std_vector_electron_sigmaIetaIeta  -> clear();  
   std_vector_electron_HoE            -> clear();            
   std_vector_electron_d0             -> clear();             
   std_vector_electron_dz             -> clear();             
   std_vector_electron_ooEooP         -> clear();         
   std_vector_electron_passConversion -> clear(); 
   std_vector_lepton_chargedHadronIso -> clear(); 
   std_vector_lepton_photonIso        -> clear();        
   std_vector_lepton_neutralHadronIso -> clear(); 
   std_vector_lepton_sumPUPt          -> clear();          
   std_vector_lepton_chargedHadronIso -> clear(); 
   std_vector_lepton_photonIso        -> clear();        
   std_vector_lepton_neutralHadronIso -> clear(); 
   std_vector_electron_effectiveArea  -> clear();  
   
} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//



void WWAnalysisSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master

  tree = FindOutput<TTree*>("nt");
  ///*** 1D histos ***/// 

  h_N_PV  = FindOutput<TH1F*>("h_N_PV");  
  
  // Counting histograms
  //----------------------------------------------------------------------------
  
  hWTrigger     = FindOutput<TH1F*>("hWTrigger");
  hWMetCut      = FindOutput<TH1F*>("hWMetCut");
  hWLowMinv     = FindOutput<TH1F*>("hWLowMinv");
  hWZVeto       = FindOutput<TH1F*>("hWZVeto");
  hWpMetCut     = FindOutput<TH1F*>("hWpMetCut");
  hWJetVeto     = FindOutput<TH1F*>("hWJetVeto");
  hWnJets       = FindOutput<TH1F*>("hWnJets");
  hWeffnJets    = FindOutput<TH1F*>("hWeffnJets");

  hWnBtaggedJets         = FindOutput<TH1F*>("hWnBtaggedJets");
  hWeffnBtaggedJets      = FindOutput<TH1F*>("hWeffnBtaggedJets");
  hWnJetsBveto           = FindOutput<TH1F*>("hWnJetsBveto");
  hWeffnJetsBveto        = FindOutput<TH1F*>("hWeffnJetsBveto");
  hNjetsTwoLeptonsLevel  = FindOutput<TH1F*>("hNjetsTwoLeptonsLevel");
  hNjetsPlot1TwoLeptonsLevel  = FindOutput<TH1F*>("hNjetsPlot1TwoLeptonsLevel");
  hNjetsPlot2TwoLeptonsLevel  = FindOutput<TH1F*>("hNjetsPlot2TwoLeptonsLevel");
  hWnJetsBvetoAfterHt    = FindOutput<TH1F*>("hWnJetsBvetoAfterHt");

  hWDeltaPhiJet = FindOutput<TH1F*>("hWDeltaPhiJet");
  hWSoftMuVeto  = FindOutput<TH1F*>("hWSoftMuVeto");
  hWExtraLepton = FindOutput<TH1F*>("hWExtraLepton");
  hWPtll        = FindOutput<TH1F*>("hWPtll");
  hWTopTagging  = FindOutput<TH1F*>("hWTopTagging");

  hWeffTrigger     = FindOutput<TH1F*>("hWeffTrigger");
  hWeffMetCut      = FindOutput<TH1F*>("hWeffMetCut");
  hWeffLowMinv     = FindOutput<TH1F*>("hWeffLowMinv");
  hWeffZVeto       = FindOutput<TH1F*>("hWeffZVeto");
  hWeffpMetCut     = FindOutput<TH1F*>("hWeffpMetCut");
  hWeffJetVeto     = FindOutput<TH1F*>("hWeffJetVeto");
  hWeffDeltaPhiJet = FindOutput<TH1F*>("hWeffDeltaPhiJet");
  hWeffSoftMuVeto  = FindOutput<TH1F*>("hWeffSoftMuVeto");
  hWeffExtraLepton = FindOutput<TH1F*>("hWeffExtraLepton");
  hWeffPtll        = FindOutput<TH1F*>("hWeffPtll");
  hWeffTopTagging  = FindOutput<TH1F*>("hWeffTopTagging");

  hLooseIso = FindOutput<TH1F*>("hLooseIso");

  // WW level histograms
  //----------------------------------------------------------------------------                                                               
  char name[80];

  for (Int_t qq = 0; qq < 4; ++qq){
    sprintf(name,"hPtLepton1WWLevel%.1i",qq);
    hPtLepton1WWLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtLepton2WWLevel%.1i",qq);
    hPtLepton2WWLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtDiLeptonWWLevel%.1i",qq);
    hPtDiLeptonWWLevel[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hMinvWWLevel%.1i",qq);
    hMinvWWLevel[qq]            = FindOutput<TH1F*>(name);
    sprintf(name,"hMtWWLevel%.1i",qq);
    hMtWWLevel[qq]              = FindOutput<TH1F*>(name);
    sprintf(name,"hpfMetWWLevel%.1i",qq);
    hpfMetWWLevel[qq]           = FindOutput<TH1F*>(name);
    sprintf(name,"hpminMetWWLevel%.1i",qq);
    hpminMetWWLevel[qq]         = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaRLeptonsWWLevel%.1i",qq);
    hDeltaRLeptonsWWLevel[qq]   = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaPhiLeptonsWWLevel%.1i",qq);
    hDeltaPhiLeptonsWWLevel[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hDPhiPtllJetWWLevel%.1i",qq);
    hDPhiPtllJetWWLevel[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hSigEl%.1i",qq);
    hSigEl[qq]                  = FindOutput<TH1F*>(name);
    sprintf(name,"hSigMu%.1i",qq);
    hSigMu[qq]                  = FindOutput<TH1F*>(name);

    sprintf(name,"hPtLepton1WWLevelNoHt%.1i",qq);
    hPtLepton1WWLevelNoHt[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtLepton2WWLevelNoHt%.1i",qq);
    hPtLepton2WWLevelNoHt[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtDiLeptonWWLevelNoHt%.1i",qq);
    hPtDiLeptonWWLevelNoHt[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hMinvWWLevelNoHt%.1i",qq);
    hMinvWWLevelNoHt[qq]            = FindOutput<TH1F*>(name);
    sprintf(name,"hMtWWLevelNoHt%.1i",qq);
    hMtWWLevelNoHt[qq]              = FindOutput<TH1F*>(name);
    sprintf(name,"hpfMetWWLevelNoHt%.1i",qq);
    hpfMetWWLevelNoHt[qq]           = FindOutput<TH1F*>(name);
    sprintf(name,"hpminMetWWLevelNoHt%.1i",qq);
    hpminMetWWLevelNoHt[qq]         = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaRLeptonsWWLevelNoHt%.1i",qq);
    hDeltaRLeptonsWWLevelNoHt[qq]   = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaPhiLeptonsWWLevelNoHt%.1i",qq);
    hDeltaPhiLeptonsWWLevelNoHt[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hDPhiPtll_JetWWLevelNoHt%.1i",qq);
    hDPhiPtllJetWWLevelNoHt[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hSigEl%.1i",qq);
    hSigElNoHt[qq]                  = FindOutput<TH1F*>(name);
    sprintf(name,"hSigMu%.1i",qq);
    hSigMuNoHt[qq]                  = FindOutput<TH1F*>(name);

    sprintf(name,"hPtLepton1WWLevelHtPlus%.1i",qq);
    hPtLepton1WWLevelHtPlus[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtLepton2WWLevelHtPlus%.1i",qq);
    hPtLepton2WWLevelHtPlus[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtDiLeptonWWLevelHtPlus%.1i",qq);
    hPtDiLeptonWWLevelHtPlus[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hMinvWWLevelHtPlus%.1i",qq);
    hMinvWWLevelHtPlus[qq]            = FindOutput<TH1F*>(name);
    sprintf(name,"hMtWWLevelHtPlus%.1i",qq);
    hMtWWLevelHtPlus[qq]              = FindOutput<TH1F*>(name);
    sprintf(name,"hpfMetWWLevelHtPlus%.1i",qq);
    hpfMetWWLevelHtPlus[qq]           = FindOutput<TH1F*>(name);
    sprintf(name,"hpminMetWWLevelHtPlus%.1i",qq);
    hpminMetWWLevelHtPlus[qq]         = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaRLeptonsWWLevelHtPlus%.1i",qq);
    hDeltaRLeptonsWWLevelHtPlus[qq]   = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaPhiLeptonsWWLevelHtPlus%.1i",qq);
    hDeltaPhiLeptonsWWLevelHtPlus[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hDPhiPtll_JetWWLevelHtPlus%.1i",qq);
    hDPhiPtllJetWWLevelHtPlus[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hSigEl%.1i",qq);
    hSigElHtPlus[qq]                  = FindOutput<TH1F*>(name);
    sprintf(name,"hSigMu%.1i",qq);
    hSigMuHtPlus[qq]                  = FindOutput<TH1F*>(name);

    sprintf(name,"hHt%.1i",qq);
    hHt[qq]                     = FindOutput<TH1F*>(name);
    sprintf(name,"hHtAfter%.1i",qq);
    hHtAfter[qq]                = FindOutput<TH1F*>(name);
  }

  // TwoLeptons level histograms
  //---------------------------------------------------------------------------

  hPtLepton1TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton1TwoLeptonsLevel");
  hPtLepton2TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton2TwoLeptonsLevel");
  hPtDiLeptonTwoLeptonsLevel      = FindOutput<TH1F*>("hPtDiLeptonTwoLeptonsLevel");
  hMinvTwoLeptonsLevel            = FindOutput<TH1F*>("hMinvTwoLeptonsLevel");
  hMtTwoLeptonsLevel              = FindOutput<TH1F*>("hMtTwoLeptonsLevel");
  hpfMetTwoLeptonsLevel           = FindOutput<TH1F*>("hpfMetTwoLeptonsLevel");
  hpminMetTwoLeptonsLevel         = FindOutput<TH1F*>("hpminMetTwoLeptonsLevel");
  hDeltaRLeptonsTwoLeptonsLevel   = FindOutput<TH1F*>("hDeltaRLeptonsTwoLeptonsLevel");
  hDeltaPhiLeptonsTwoLeptonsLevel = FindOutput<TH1F*>("hDeltaPhiLeptonsTwoLeptonsLevel");
  hDPhiPtllJetTwoLeptonsLevel     = FindOutput<TH1F*>("hDPhiPtllJetTwoLeptonsLevel");
  hSigMuNoHtTwoLeptonsLevel       = FindOutput<TH1F*>("hSigMuNoHtTwoLeptonsLevel");
  hSigElNoHtTwoLeptonsLevel       = FindOutput<TH1F*>("hSigElNoHtTwoLeptonsLevel");
}

//------------------------------------------------------------------------------
// IsTightLepton
//
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
// egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium
//------------------------------------------------------------------------------
bool WWAnalysisSelector::IsTightLepton(int k)
{
  bool is_tight_lepton = false;

  // Muon ID
  if (fabs(std_vector_lepton_id->at(k)) == 13){

    if (_MuonID.Contains("MediumID")){
      if (std_vector_lepton_isMediumMuon->at(k) == 1){
	is_tight_lepton = true;
	if(_MuonID == "MediumIDTighterIP"){
	  is_tight_lepton = false;
	  if (fabs(std_vector_lepton_BestTrackdxy -> at(k)) < 0.02)
	    if (fabs(std_vector_lepton_BestTrackdz -> at(k)) < 0.1)
	      is_tight_lepton = true;
	}
      }
    }
    
    else if (_MuonID.Contains("TightID") ){
      if (std_vector_lepton_isTightMuon->at(k) == 1){
	is_tight_lepton = true;
	if(_MuonID == "TightIDTighterIP"){
	  is_tight_lepton = false;
	  if (fabs(std_vector_lepton_BestTrackdxy -> at(k)) < 0.02)
	    if (fabs(std_vector_lepton_BestTrackdz -> at(k)) < 0.1)
	      is_tight_lepton = true;
	}
      }
    }

  }
  
  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_id->at(k)) == 11)
    {
      float aeta = fabs(std_vector_electron_scEta->at(k));
      
      if (aeta <= 1.479)
	{
	  if (fabs(std_vector_electron_deltaEtaIn->at(k)) < 0.008925 &&
	      fabs(std_vector_electron_deltaPhiIn->at(k)) < 0.035973 &&
	      std_vector_electron_sigmaIetaIeta->at(k)    < 0.009996 &&
	      std_vector_electron_HoE->at(k)              < 0.050537 &&
	      fabs(std_vector_electron_d0->at(k))         < 0.012235 &&
	      fabs(std_vector_electron_dz->at(k))         < 0.042020 &&
	      fabs(std_vector_electron_ooEooP->at(k))     < 0.091942 &&
	      ElectronIsolation(k)                        < 0.107587 &&
	      !std_vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
      else if (aeta > 1.479 && aeta < 2.5)
	{
	  if (fabs(std_vector_electron_deltaEtaIn->at(k)) < 0.007429 &&
	      fabs(std_vector_electron_deltaPhiIn->at(k)) < 0.067879 &&
	      std_vector_electron_sigmaIetaIeta->at(k)    < 0.030135 &&
	      std_vector_electron_HoE->at(k)              < 0.086782 &&
	      fabs(std_vector_electron_d0->at(k))         < 0.036719 &&
	      fabs(std_vector_electron_dz->at(k))         < 0.138142 &&
	      fabs(std_vector_electron_ooEooP->at(k))     < 0.100683 &&
	      ElectronIsolation(k)                        < 0.113254 &&
	      !std_vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
    }

  return is_tight_lepton;

}


//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float WWAnalysisSelector::MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_id->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 13) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      0.5*std_vector_lepton_sumPUPt->at(k)));

  relative_isolation /= pt;

  return relative_isolation;
}


//------------------------------------------------------------------------------
// ElectronIsolation
//------------------------------------------------------------------------------
float WWAnalysisSelector::ElectronIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_id->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 11) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      jetRho*std_vector_electron_effectiveArea->at(k)));
  
  relative_isolation /= pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// IsIsolatedLepton
//------------------------------------------------------------------------------
bool WWAnalysisSelector::IsIsolatedLepton(int k)
{
  float id = std_vector_lepton_id->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}
