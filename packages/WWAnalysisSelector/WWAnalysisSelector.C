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

  //Original 8 TeV DY-MVA Variables
  tree->Branch("fullpmet",&fullpmet,"fullpmet");
  tree->Branch("trkpmet",&trkpmet,"trkpmet");
  tree->Branch("ratioMet",&ratioMet,"ratioMet");
  tree->Branch("ptll",&ptll,"ptll");
  tree->Branch("mth",&mth,"mth");
  tree->Branch("jetpt1",&jetpt1,"jetpt1");
  tree->Branch("ptWW",&ptWW,"ptWW");
  tree->Branch("dphilljet",&dphilljet,"dphilljet");
  tree->Branch("dphillmet",&dphillmet,"dphillmet");
  tree->Branch("dphijet1met",&dphijet1met,"dphijet1met");
  tree->Branch("nvtx",&nvtx,"nvtx");

  //Variables Added to Be Studied
  tree->Branch("pt1",&pt1,"pt1");
  tree->Branch("pt2",&pt2,"pt2");
  tree->Branch("ptll",&ptll,"ptll");
  tree->Branch("drll",&drll,"drll");
  tree->Branch("pfType1Met",&pfType1Met,"pfType1Met");
  tree->Branch("trkMet",&trkMet,"trkMet");
  tree->Branch("dphill",&dphill,"dphill");
  tree->Branch("jetphi1",&jetphi1,"jetphi1");
  tree->Branch("pfType1Metphi",&pfType1Metphi,"pfType1Metphi");
  tree->Branch("mll",&mll,"mll");
  tree->Branch("mpmet",&mpmet,"mpmet");
  tree->Branch("metvar",&metvar,"metvar");
  tree->Branch("njet",&njet,"njet");
  tree->Branch("dphilmet1",&dphilmet1,"dphilmet1");
  tree->Branch("dphilmet2",&dphilmet2,"dphilmet2");
  tree->Branch("bveto_ip",&bveto_ip,"bveto_ip");
  tree->Branch("nbjettche",&nbjettche,"nbjettche");
  tree->Branch("Ht",&Ht,"Ht");
  tree->Branch("baseW",&baseW,"baseW");

  //------------------------------------------------------------------------------
  // Create histos
  //------------------------------------------------------------------------------
  
  h_n_PV = CreateH1F("h_n_PV","h_n_PV",50,0.,10.);

  //Z + Jets Data Driven histograms - nMetCut x nJet
  //----------------------------------------------------------------------------       
  
  for (unsigned int j = 0; j < 4; ++j){
    for (size_t nC=0; nC<numberMetCuts; nC++) {
      hNinZevents     [nC][j] = CreateH1F(Form("hNinZevents%.1i%.1i",      MetCut[nC],j), "",    3, 0,    3);
      hNoutZevents    [nC][j] = CreateH1F(Form("hNoutZevents%.1i%.1i",     MetCut[nC],j), "",    3, 0,    3);
      hNinLooseZevents[nC][j] = CreateH1F(Form("hNinLooseZevents%.1i%.1i", MetCut[nC],j), "",    3, 0,    3);
      hMassInZevents  [nC][j] = CreateH1F(Form("hMassInZevents%.1i%.1i",   MetCut[nC],j), "", 3000, 0, 3000);
      hMassOutZevents [nC][j] = CreateH1F(Form("hMassOutZevents%.1i%.1i",  MetCut[nC],j), "", 3000, 0, 3000);
    } 
  }

  // Counting histograms    
  //---------------------------------------------------------------------------- 
  
  for (Int_t nC=0; nC<4; nC++) {

    hWTrigger[nC]     = CreateH1F(Form("hWTrigger%.1i", nC),                        "", 10, 0, 10);
    hWMetCut[nC]      = CreateH1F(Form("hWMetCut%.1i", nC),                         "", 10, 0, 10);
    hWLowMinv[nC]     = CreateH1F(Form("hWLowMinv%.1i", nC),                        "", 10, 0, 10);
    hWZVeto[nC]       = CreateH1F(Form("hWZVeto%.1i", nC),                          "", 10, 0, 10);
    hWpMetCut[nC]     = CreateH1F(Form("hWpMetCut%.1i", nC),                        "", 10, 0, 10);
    hWJetVeto[nC]     = CreateH1F(Form("hWJetVeto%.1i", nC),                        "", 10, 0, 10);
    hWnJets[nC]       = CreateH1F(Form("hWnJets%.1i", nC),                          "", 10, 0, 10);
    hWeffnJets[nC]    = CreateH1F(Form("hWeffnJets%.1i", nC),                       "", 10, 0, 10);
    hWnBtaggedJets[nC]     = CreateH1F(Form("hWnBtaggedJets%.1i", nC),              "", 10, 0, 10);
    hWeffnBtaggedJets[nC] = CreateH1F(Form("hWeffnBtaggedJets%.1i", nC),            "", 10, 0, 10);
    hWnJetsBveto[nC]    = CreateH1F(Form("hWnJetsBveto%.1i", nC),                   "", 10, 0, 10);
    hWeffnJetsBveto[nC] = CreateH1F(Form("hWeffnJetsBveto%.1i", nC),                "", 10, 0, 10);
    hNjetsTwoLeptonsLevel[nC]    = CreateH1F(Form("hNjetsTwoLeptonsLevel%.1i", nC), "", 10, 0, 10);
    hWeffnJetsBvetoAfterHt[nC] = CreateH1F(Form("hWeffnJetsBvetoAfterHt%.1i", nC),  "", 10, 0, 10);
    
    hWDeltaPhiJet[nC] = CreateH1F(Form("hWDeltaPhiJet%.1i", nC), "", 10, 0, 10);
    hWSoftMuVeto[nC]  = CreateH1F(Form("hWSoftMuVeto%.1i", nC),  "", 10, 0, 10);
    hWExtraLepton[nC] = CreateH1F(Form("hWExtraLepton%.1i", nC), "", 10, 0, 10);
    hWPtll[nC]        = CreateH1F(Form("hWPtll%.1i", nC),        "", 10, 0, 10);
    hWTopTagging[nC]  = CreateH1F(Form("hWTopTagging%.1i", nC),  "", 10, 0, 10);
    
    hWeffTrigger[nC]     = CreateH1F(Form("hWeffTrigger%.1i", nC),     "", 10, 0, 10);
    hWeffMetCut[nC]      = CreateH1F(Form("hWeffMetCut%.1i", nC),      "", 10, 0, 10);
    hWeffLowMinv[nC]     = CreateH1F(Form("hWeffLowMinv%.1i", nC),     "", 10, 0, 10);
    hWeffZVeto[nC]       = CreateH1F(Form("hWeffZVeto%.1i", nC),       "", 10, 0, 10);
    hWeffpMetCut[nC]     = CreateH1F(Form("hWeffpMetCut%.1i", nC),     "", 10, 0, 10);
    hWeffJetVeto[nC]     = CreateH1F(Form("hWeffJetVeto%.1i", nC),     "", 10, 0, 10);
    hWeffDeltaPhiJet[nC] = CreateH1F(Form("hWeffDeltaPhiJet%.1i", nC), "", 10, 0, 10);
    hWeffSoftMuVeto[nC]  = CreateH1F(Form("hWeffSoftMuVeto%.1i", nC),  "", 10, 0, 10);
    hWeffExtraLepton[nC] = CreateH1F(Form("hWeffExtraLepton%.1i", nC), "", 10, 0, 10);
    hWeffPtll[nC]        = CreateH1F(Form("hWeffPtll%.1i", nC),        "", 10, 0, 10);
    hWeffTopTagging[nC]  = CreateH1F(Form("hWeffTopTagging%.1i", nC),  "", 10, 0, 10);
    
    hLooseIso[nC] = CreateH1F(Form("hLooseIso%.1i", nC),  "", 100, 0, 10);
    
    // TwoLeptons level histograms   
    //----------------------------------------------------------------------------
    
    hNVtx[nC]                           = CreateH1F(Form("hNVtx%.1i", nC),                           "",  100,0., 100);  
    hPtLepton1TwoLeptonsLevel[nC]       = CreateH1F(Form("hPtLepton1TwoLeptonsLevel%.1i", nC),       "", 3000,0.,3000);
    hPtLepton2TwoLeptonsLevel[nC]       = CreateH1F(Form("hPtLepton2TwoLeptonsLevel%.1i", nC),       "", 3000,0.,3000);
    hPtDiLeptonTwoLeptonsLevel[nC]      = CreateH1F(Form("hPtDiLeptonTwoLeptonsLevel%.1i", nC),      "", 3000,0.,3000);
    hMinvTwoLeptonsLevel[nC]            = CreateH1F(Form("hMinvTwoLeptonsLevel%.1i", nC),            "", 3000,0.,3000);
    hMtTwoLeptonsLevel[nC]              = CreateH1F(Form("hMtTwoLeptonsLevel%.1i", nC),              "", 3000,0.,3000);
    hpfMetTwoLeptonsLevel[nC]           = CreateH1F(Form("hpfMetTwoLeptonsLevel%.1i", nC),           "", 3000,0.,3000);
    htrkMetTwoLeptonsLevel[nC]          = CreateH1F(Form("htrkMetTwoLeptonsLevel%.1i", nC),          "", 3000,0.,3000);
    hpminMetTwoLeptonsLevel[nC]         = CreateH1F(Form("hpminMetTwoLeptonsLevel%.1i", nC),         "", 3000,0.,3000);
    hDeltaRLeptonsTwoLeptonsLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsTwoLeptonsLevel%.1i", nC),   "",  50, 0.,   5);
    hDeltaPhiLeptonsTwoLeptonsLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsTwoLeptonsLevel%.1i", nC), "",  32, 0., 3.2);
    hDPhiPtllJetTwoLeptonsLevel[nC]     = CreateH1F(Form("hDPhiPtllJetTwoLeptonsLevel%.1i", nC),     "",  32, 0., 3.2);
    hNjetsPlot1TwoLeptonsLevel[nC]      = CreateH1F(Form("hNjetsPlot1TwoLeptonsLevel%.1i", nC),      "",  10, 0.,  10);
    hNjetsPlot2TwoLeptonsLevel[nC]      = CreateH1F(Form("hNjetsPlot2TwoLeptonsLevel%.1i", nC),      "",  10, 0.,  10);
    hSigMuNoHtTwoLeptonsLevel[nC]       = CreateH1F(Form("hSigMuNoHtTwoLeptonsLevel%.1i", nC),       "",  10, 0.,  10);
    hSigElNoHtTwoLeptonsLevel[nC]       = CreateH1F(Form("hSigElNoHtTwoLeptonsLevel%.1i", nC),       "",  10, 0.,  10);
    hDxyTwoLeptonsLevel[nC]             = CreateH1F(Form("hDxyTwoLeptonsLevel%.1i", nC),             "", 1000,-0.1,0.1);
    hDzTwoLeptonsLevel[nC]              = CreateH1F(Form("hDzTwoLeptonsLevel%.1i", nC),              "", 1000,-0.5,0.5);
    
    hsoftMuPt[nC]                       = CreateH1F(Form("hsoftMuPt%.1i", nC),                       "", 3000,0.,3000);
    hjetPt[nC]                          = CreateH1F(Form("hjetPt%.1i", nC),                          "", 3000,0.,3000);
    
    // WW level histograms     
    //----------------------------------------------------------------------------  
    
    hWnJetsBvetoAfterHt[nC] = CreateH1F(Form("hWnJetsBvetoAfterHt%.1i", nC),  "", 10, 0, 10); 
    
    hPtLepton1WWLevel[nC]       = CreateH1F(Form("hPtLepton1WWLevel%.1i", nC),        "", 3000,0.,3000);
    hPtLepton2WWLevel[nC]       = CreateH1F(Form("hPtLepton2WWLevel%.1i", nC),        "", 3000,0.,3000);
    hPtDiLeptonWWLevel[nC]      = CreateH1F(Form("hPtDiLeptonWWLevel%.1i", nC),       "", 3000,0.,3000);
    hMinvWWLevel[nC]            = CreateH1F(Form("hMinvWWLevel%.1i", nC),             "", 3000,0.,3000);
    hMtWWLevel[nC]              = CreateH1F(Form("hMtWWLevel%.1i", nC),               "", 3000,0.,3000);
    hpfMetWWLevel[nC]           = CreateH1F(Form("hpfMetWWLevel%.1i", nC),            "", 3000,0.,3000);
    htrkMetWWLevel[nC]          = CreateH1F(Form("htrkMetWWLevel%.1i", nC),           "", 3000,0.,3000);
    hpminMetWWLevel[nC]         = CreateH1F(Form("hpminMetWWLevel%.1i", nC),          "", 3000,0.,3000);
    hDeltaRLeptonsWWLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),    "",   50, 0,   5);
    hDeltaPhiLeptonsWWLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC),  "",   32, 0, 3.2);
    hDPhiPtllJetWWLevel[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),      "",   32, 0, 3.2);
    hSigEl[nC]                  = CreateH1F(Form("hSigEl%.1i", nC),                   "", 3000,0.,3000);
    hSigMu[nC]                  = CreateH1F(Form("hSigMu%.1i", nC),                   "", 3000,0.,3000);

    hHt[nC]                     = CreateH1F(Form("hHt%.1i",               nC),        "", 3000, 0, 3000);
    hHtAfter[nC]                = CreateH1F(Form("hHtAfter%.1i",          nC),        "", 3000, 0, 3000);
    
    // BVeto level histograms 
    //----------------------------------------------------------------------------
    
    hbVetoOld[nC]                       = CreateH1F(Form("hbVetoOld%.1i", nC),                       "", 3000, 0, 3000);
    hbVetoMu[nC]                        = CreateH1F(Form("hbVetoMu%.1i", nC),                        "", 3000, 0, 3000);
    hbVetoCsvv2ivfLoose[nC]             = CreateH1F(Form("hbVetoCsvv2ivfLoose%.1i", nC),             "", 3000, 0, 3000);
    hbVetoCsvv2ivfMedium[nC]            = CreateH1F(Form("hbVetoCsvv2ivfMedium%.1i", nC),            "", 3000, 0, 3000);
    hbVetoCsvv2ivfTight[nC]             = CreateH1F(Form("hbVetoCsvv2ivfTight%.1i", nC),             "", 3000, 0, 3000);
    hbVetoCsvv2ivfLooseAndMu[nC]        = CreateH1F(Form("hbVetoCsvv2ivfLooseAndMu%.1i", nC),        "", 3000, 0, 3000);
    hbVetoCsvv2ivfRecommended[nC]       = CreateH1F(Form("hbVetoCsvv2ivfRecommended%.1i", nC),       "", 3000, 0, 3000);
    hbVetoCsvv2ivfRecommendedAndMu[nC]  = CreateH1F(Form("hbVetoCsvv2ivfRecommendedAndMu%.1i", nC),  "", 3000, 0, 3000);
    
    //Control Region Plots
    //----------------------------------------------------------------------------                
    
    hpfMetCR[nC]  = CreateH1F(Form("hpfMetCR%.1i", nC),          "", 3000, 0, 3000);
    htrkMetCR[nC] = CreateH1F(Form("htrkMetCR%.1i", nC),         "", 3000, 0, 3000);
    hmpMetCR[nC]  = CreateH1F(Form("hmpMetCR%.1i", nC),          "", 3000, 0, 3000);
    
    // monoH level histograms     
    //----------------------------------------------------------------------------  
    
    for (Int_t jj=0; jj<4; jj++) {
      
      hPtLepton1monoHLevel[nC][jj]       = CreateH1F(Form("hPtLepton1monoHLevel%.1i%.1i", nC,jj),       "", 3000, 0, 3000);
      hPtLepton2monoHLevel[nC][jj]       = CreateH1F(Form("hPtLepton2monoHLevel%.1i%.1i", nC,jj),       "", 3000, 0, 3000);
      hPtDiLeptonmonoHLevel[nC][jj]      = CreateH1F(Form("hPtDiLeptonmonoHLevel%.1i%.1i", nC,jj),      "", 3000, 0, 3000);
      hMinvmonoHLevel[nC][jj]            = CreateH1F(Form("hMinvmonoHLevel%.1i%.1i", nC,jj),            "", 3000, 0, 3000);
      hMtmonoHLevel[nC][jj]              = CreateH1F(Form("hMtmonoHLevel%.1i%.1i", nC,jj),              "", 3000, 0, 3000);
      hMt1monoHLevel[nC][jj]             = CreateH1F(Form("hMt1monoHLevel%.1i%.1i", nC,jj),             "", 3000, 0, 3000);
      hMt2monoHLevel[nC][jj]             = CreateH1F(Form("hMt2monoHLevel%.1i%.1i", nC,jj),             "", 3000, 0, 3000);
      hpfMetmonoHLevel[nC][jj]           = CreateH1F(Form("hpfMetmonoHLevel%.1i%.1i", nC,jj),           "", 3000, 0, 3000);
      hpminMetmonoHLevel[nC][jj]         = CreateH1F(Form("hpminMetmonoHLevel%.1i%.1i", nC,jj),         "", 3000, 0, 3000);
      hDeltaRLeptonsmonoHLevel[nC][jj]   = CreateH1F(Form("hDeltaRLeptonsmonoHLevel%.1i%.1i", nC,jj),   "",   50, 0,    5);
      hDeltaPhiLeptonsmonoHLevel[nC][jj] = CreateH1F(Form("hDeltaPhiLeptonsmonoHLevel%.1i%.1i", nC,jj), "",   32, 0,  3.2);
      hDPhiPtllJetmonoHLevel[nC][jj]     = CreateH1F(Form("hDPhiPtllJetmonoHLevel%.1i%.1i", nC,jj),     "",   32, 0,  3.2);
      hDPhillMetmonoHLevel[nC][jj]       = CreateH1F(Form("hDPhillMetmonoHLevel%.1i%.1i", nC,jj),       "",   32, 0,  3.2);
      hPtWWmonoHLevel[nC][jj]            = CreateH1F(Form("hPtWWmonoHLevel%.1i%.1i", nC,jj),            "", 3000, 0, 3000);
      hMcmonoHLevel[nC][jj]              = CreateH1F(Form("hMcmonoHLevel%.1i%.1i", nC,jj),              "", 3000, 0, 3000);
      hTrkMetmonoHLevel[nC][jj]          = CreateH1F(Form("hTrkMetmonoHLevel%.1i%.1i", nC,jj),          "", 3000, 0, 3000);
      hHtmonoHLevel[nC][jj]              = CreateH1F(Form("hHtmonoHLevel%.1i%.1i", nC,jj),              "", 3000, 0, 3000);
      hdphijetjetmonoHLevel[nC][jj]      = CreateH1F(Form("hdphijetjetmonoHLevel%.1i%.1i", nC,jj),      "",   32, 0,  3.2);
      hdetajetjetmonoHLevel[nC][jj]      = CreateH1F(Form("hdetajetjetmonoHLevel%.1i%.1i", nC,jj),      "",   50, 0,  5.0);
    }
  }
}

void WWAnalysisSelector::InsideLoop() {
  
  //Assigning values to vectors
  //----------------------------------------------------------------------------
  Assign("std_vector_lepton_pt",std_vector_lepton_pt);
  Assign("std_vector_jet_pt",std_vector_jet_pt);    
  Assign("std_vector_lepton_muSIP3D",std_vector_lepton_muSIP3D);
  Assign("std_vector_lepton_elSIP3D",std_vector_lepton_elSIP3D);
  Assign("std_vector_lepton_flavour",std_vector_lepton_flavour);     
  Assign("std_vector_lepton_isTightMuon",std_vector_lepton_isTightMuon);
  Assign("std_vector_lepton_eleIdMedium",std_vector_lepton_eleIdMedium);
  Assign("std_vector_lepton_eleIdTight",std_vector_lepton_eleIdTight);
  Assign("std_vector_lepton_chargedHadronIso",std_vector_lepton_chargedHadronIso);
  Assign("std_vector_lepton_photonIso",std_vector_lepton_photonIso);       
  Assign("std_vector_lepton_neutralHadronIso",std_vector_lepton_neutralHadronIso);
  Assign("std_vector_lepton_sumPUPt",std_vector_lepton_sumPUPt); 
  Assign("std_vector_electron_effectiveArea",std_vector_electron_effectiveArea); 
  Assign("std_vector_lepton_BestTrackdxy",std_vector_lepton_BestTrackdxy);
  Assign("std_vector_lepton_BestTrackdz",std_vector_lepton_BestTrackdz);  
  Assign("std_vector_lepton_phi",std_vector_lepton_phi);
  Assign("std_vector_lepton_eta",std_vector_lepton_eta);
  Assign("std_vector_jet_eta",std_vector_jet_eta);
  Assign("std_vector_jet_phi",std_vector_jet_phi);  
  Assign("std_vector_lepton_isMediumMuon",std_vector_lepton_isMediumMuon);
  Assign("std_vector_jet_csvv2ivf",std_vector_jet_csvv2ivf);
  Assign("std_vector_jet_softMuPt",std_vector_jet_softMuPt);
  Assign("std_vector_electron_hOverE",std_vector_electron_hOverE);
  Assign("std_vector_electron_dEtaIn",std_vector_electron_dEtaIn);
  Assign("std_vector_electron_dPhiIn",std_vector_electron_dPhiIn);
  Assign("std_vector_electron_ooEmooP",std_vector_electron_ooEmooP);

  //Assigning values to float variables
  //----------------------------------------------------------------------------
  if (_Signal != "VV50" && _Signal != "Data201550" && _Signal != "Top50")
    Assign("GEN_weight_SM",GEN_weight_SM);
  else 
    GEN_weight_SM = 1;
  if (_Signal == "WW50")
    Assign("nllW",nllW);
  Assign("phi1",phi1);
  Assign("phi2",phi2);
  Assign("puW",puW);  
  Assign("effW",effW);
  Assign("triggW",triggW);
  Assign("trigger",trigger);
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
  Assign("jeteta1",jeteta1);
  Assign("jetpt2",jetpt2);
  Assign("jetphi2",jetphi2);
  Assign("jeteta2",jeteta2);
  Assign("pfType1Metphi",pfType1Metphi);
  Assign("dphillmet",dphillmet);
  Assign("nvtx",nvtx);

  //Assigning values to integer variables
  //----------------------------------------------------------------------------
  Assign("bveto_ip",bveto_ip);

  //Creating the variables we need
  //----------------------------------------------------------------------------
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  //self-definition of njet!!!
  //njet = 0;
  //for (unsigned int i = 0; i < std_vector_lepton_pt.size(); ++i)
  //if (std_vector_lepton_pt.at(i) > 30)
  //++njet;

  //NNLL corrections
  if (_Signal == "WW50")
    baseW = baseW * nllW * (12.178 / 10.481);

  //NNLO corrections
  if (_Signal == "WW25")
    baseW = baseW * (12.178 / 10.481);

  //Correcting for the Negative Weight Scaling  
  if (_Signal == "WJets50")
    baseW = ( (GEN_weight_SM/abs(GEN_weight_SM)) / 0.68394 ) * baseW;

  else if (_Signal == "DY50"){
    baseW = ( (GEN_weight_SM/abs(GEN_weight_SM)) / 0.670032 ) * baseW;
    //if (channel == 1) //EE
    //baseW = baseW * 1.07921;
    //else if (channel == 0) //MuMu
    //baseW = baseW * 0.85099;
  }

  else if (_Signal == "SingleTop50")
    baseW = ( (GEN_weight_SM/abs(GEN_weight_SM)) / 0.215131 ) * baseW;

  else if (_Signal == "TTJets50")
    baseW = ( (GEN_weight_SM/abs(GEN_weight_SM)) / 0.331907 ) * baseW;

  else if (_Signal == "DY25"){
    baseW = ( (GEN_weight_SM/abs(GEN_weight_SM)) / 0.72760 ) * baseW;
    //if (channel == 1) //EE
    //baseW = baseW * 1.07921;
    //else if (channel == 0) //MuMu
    //baseW = baseW * 0.85099;
  }

  else
    baseW = baseW * 1.;
 
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * _Luminosity;
  
  if (_Signal.Contains("Data"))
    totalW = 1.0;

  h_n_PV -> Fill(1,efficiencyW);  
  
  dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  int jetbin = njet;
  
  //Building mpmet
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  
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
  
  //defining nextra
  nextra = 0;
  for (unsigned int i = 0; i < std_vector_lepton_pt.size(); ++i)
    if(std_vector_lepton_pt.at(i) > 10)
      ++nextra;
  nextra = nextra - 2;
       
  Float_t dRjet1  = 100.;
  Float_t dRjet2  = 100.;

  //building Mt of the two Ws
  if( pfType1Met > 0){
    if( pt1 > 0 )
      Mt1 = sqrt(2*pt1*pfType1Met*(1 - cos(dphilmet1)));  
    if( pt2 > 0 )
      Mt2 = sqrt(2*pt2*pfType1Met*(1 - cos(dphilmet2)));  
  }
  
  //building Mc
  if (ptll > 0 && mll > 0 && pfType1Met > 0)
    Mc = sqrt( pow(sqrt(ptll*ptll + mll*mll) + pfType1Met,2) - pow(ptll + pfType1Met,2) );
  
  //building dphijetjet
  if (jetpt1 > 0 && jetpt2 > 0)
    if (jetphi1 > -3.5 && jetphi1 < 3.5)
      if (jetphi2 > -3.5 && jetphi2 < 3.5){
	dphijetjet = fabs(jetphi2 - jetphi1);
	if (dphijetjet > TMath::Pi())
	  dphijetjet = 2*TMath::Pi() - dphijetjet;
      }

  //building detajetjet
  if (jetpt1 > 0 && jetpt2 > 0)
    detajetjet = fabs(jeteta1 - jeteta2);

  /*
  //building dRjet1
  Float_t dRjet1  = 100.;
  TLorentzVector vjet1(0.,0.,0.,0.);
  TLorentzVector vlep1(0.,0.,0.,0.);
  if(std_vector_lepton_pt.at(0) > 0.){
    vlep1.SetPtEtaPhiM(std_vector_lepton_pt.at(0),std_vector_lepton_eta.at(0),std_vector_lepton_phi.at(0),0.);
    for (unsigned int i = 0; i < std_vector_jet_pt.size(); ++i)
      if(std_vector_jet_pt.at(i) > 0.){
	vjet1.SetPtEtaPhiM(std_vector_jet_pt.at(i),std_vector_jet_eta.at(i),std_vector_jet_phi.at(i),0.);
	if( vlep1.DeltaR(vjet1) < dRjet1){
	  dRjet1 = vlep1.DeltaR(vjet1);
	}
      }
  }

  //building dRjet2
  Float_t dRjet2  = 100.;
  TLorentzVector vjet2(0.,0.,0.,0.);
  TLorentzVector vlep2(0.,0.,0.,0.);
  if(std_vector_lepton_pt.at(1) > 0.){
    vlep2.SetPtEtaPhiM(std_vector_lepton_pt.at(1),std_vector_lepton_phi.at(1),std_vector_lepton_eta.at(1),0.);
    for (unsigned int i = 0; i < njet; ++i)
      if(std_vector_jet_pt.at(i) > 0.){
	vjet2.SetPtEtaPhiM(std_vector_jet_pt.at(i),std_vector_jet_eta.at(i),std_vector_jet_phi.at(i),0.);
	if( vlep2.DeltaR(vjet2) < dRjet2){
	  dRjet2 = vlep2.DeltaR(vjet2);
	}
      }
  }
*/

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
  
   //Building b-veto csvv2ivf Loose (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfLoose = false;

   for (unsigned int i = 0; i < std_vector_jet_pt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 15)
       if (std_vector_jet_csvv2ivf.at(i) > 0.423)
	 bvetocsvv2ivfLoose = true;

   //Building b-veto csvv2ivf Medium (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfMedium = false;

   for (unsigned int i = 0; i < std_vector_jet_pt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 15)
       if (std_vector_jet_csvv2ivf.at(i) > 0.814)
	 bvetocsvv2ivfMedium = true;

   //Building b-veto csvv2ivf Tight (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfTight = false;

   for (unsigned int i = 0; i < std_vector_jet_pt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 15)
       if (std_vector_jet_csvv2ivf.at(i) > 0.941)
	 bvetocsvv2ivfTight = true;

   //Building b-veto csvv2ivf Recommended (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfRecommended = false;

   for (unsigned int i = 0; i < std_vector_jet_pt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 15)
       if (std_vector_jet_csvv2ivf.at(i) > 0.605)
	 bvetocsvv2ivfRecommended = true;

   //Building a Soft Muon B-Veto (Quite Rudely)
   bool bvetoMu = false;
   
   for (unsigned int i = 0; i < std_vector_jet_softMuPt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 10 && std_vector_jet_pt.at(i) < 30)
       if (std_vector_jet_softMuPt.at(i) > 3){
	 hsoftMuPt[0]->Fill(std_vector_jet_softMuPt.at(i), totalW);
	 hjetPt[0]   ->Fill(std_vector_jet_pt.at(i),       totalW);
	 bvetoMu = true;
	 break;
       }
   
   //building ptWW
   TLorentzVector L1,L2;
   TLorentzVector MET;
   ptWW = 0.;  

   if (pt1>0 && pt2>0) {  
     L1.SetPtEtaPhiM(pt1, 0, phi1, 0.);
     L2.SetPtEtaPhiM(pt2, 0, phi2, 0.);
     MET.SetPtEtaPhiM(pfType1Met, 0, pfType1Metphi, 0.);
     ptWW = (L1+L2+MET).Pt();
   }

   //Building Number of Gen Jet
   njetGen = 0.;
   /*
   for (unsigned int i = 0; i < std_vector_jetGen_pt.size(); ++i)
   if (std_vector_jetGen_pt.at(i) > 30)
   ++njetGen;       
   */
 

   //look for the two most energetics isolated and identified leptons (20GeV - 20GeV)
   unsigned int a = 100;
   for (unsigned int first = 0; first < std_vector_lepton_pt.size(); ++first)
     if (std_vector_lepton_pt.at(first) > 20)
       if (IsIsolatedLepton(first))
	 if (IsTightLepton(first,_MuonID)){
	   a = first;
	   break;
	 }

   unsigned int b = 100;
   for (unsigned int second = a + 1; second < std_vector_lepton_pt.size(); ++second)
     if (std_vector_lepton_pt.at(second) > 20)
       if (IsIsolatedLepton(second))
         if (IsTightLepton(second,_MuonID)){
	   b = second;
	   break;
	 }

//   cout<<a<<"    "<<b<<endl;


//   unsigned int a = 0;
//   unsigned int b = 1;

   //building Ht
   if(a != 100 && b != 100){
     Ht = std_vector_lepton_pt.at(a) + std_vector_lepton_pt.at(b) + pfType1Met;
   
     if(njet > 10) njet = 10;
     for (int i = 0; i < njet; ++i)
       if(std_vector_jet_pt.at(i) > 30)
	 Ht += std_vector_jet_pt.at(i);
   }

   // The selection begins here
   //--------------------------------------------------------------------------
   if (trigger == 1)
     if (std_vector_lepton_pt.at(0) > 20)
       if (std_vector_lepton_pt.at(1) > 20) 
	 if (a != 100)
	 if (b != 100)
	 if ((_SameSign == "SS" && ch1*ch2 > 0) || (_SameSign == "OS" && ch1*ch2 < 0))
	   if ( (SelectedChannel == -1)                                     || 
		(channel == SelectedChannel)                                || 
		(SelectedChannel == 4 && (channel == 2 || channel == 3) )   || 
		(SelectedChannel == 5 && (channel == 0 || channel == 1) )   ){
	     /*	    
		    if (IsTightLepton(0,_MuonID) && !IsTightLepton(1,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(1), totalW);
		    if (IsTightLepton(1,_MuonID) && !IsTightLepton(0,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(0), totalW);
	     */
	     
	     //lepton a
	     hDxyTwoLeptonsLevel[0] -> Fill(std_vector_lepton_BestTrackdxy.at(a),totalW);
	     hDzTwoLeptonsLevel[0]  -> Fill(std_vector_lepton_BestTrackdz.at(a) ,totalW);
	     
	     //lepton b
	     hDxyTwoLeptonsLevel[0] -> Fill(std_vector_lepton_BestTrackdxy.at(b),totalW);
	     hDzTwoLeptonsLevel[0]  -> Fill(std_vector_lepton_BestTrackdz.at(b) ,totalW);
	     
	     //	     if (IsIsolatedLepton(0))
	     //if (IsIsolatedLepton(1))
	     //	 if (IsTightLepton(0,_MuonID))
	     //	   if (IsTightLepton(1,_MuonID)){

	     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	     //
	     // Data Driven methods
	     //
	     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	     

	     if (pfType1Met > 20 && mpmet > 20 && mll > 12 && ptll > 45 && nextra == 0 && (dphiv || channel == 2 || channel == 3)) {
	     
	       //Z + Jets
	       //----------------------------------------------------------------------

	       if (dphiv && bveto_mu && (bveto_ip && (nbjettche == 0 || njet > 3))) {
		 
		 // Loop over the metvar bins
		 //----------------------------------------------------------------------

		 for (size_t mc=0; mc<numberMetCuts; mc ++) {
		   
		   if (metvar > MetCut[mc] && fabs(mll - ZMASS) < 7.5) {
		     hNinLooseZevents[mc][3]->Fill(1,totalW);
		     for (int jetNumber = 0; jetNumber < 3 ; ++jetNumber){
		       if (jetbin >= 3) jetbin = 2;
		       if(jetNumber == jetbin){
			 hNinLooseZevents[mc][jetNumber]->Fill(1,totalW);
		       }
		     }
		   }
		   if (metvar > MetCut[mc] && metvar < MetCut[mc+1]) {   
		     if (fabs(mll - ZMASS) < 7.5) {
		       hNinZevents[mc][3]   ->Fill(  1, totalW);
		       hMassInZevents[mc][3]->Fill(mll, totalW);
		       for (int jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			 if (jetbin >= 3) jetbin = 2;
			 if(jetNumber == jetbin){
			   hNinZevents[mc][jetNumber]   ->Fill(  1, totalW);
			   hMassInZevents[mc][jetNumber]->Fill(mll, totalW);
			 }
		       }
		     }
		     else if (fabs(mll - ZMASS) > 15) {  
		       hNoutZevents[mc][3]   ->Fill(  1, totalW);
		       hMassOutZevents[mc][3]->Fill(mll, totalW);
		       for (int jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			 if (jetbin >= 3) jetbin = 2;
			 if(jetNumber == jetbin){
			   hNoutZevents[mc][jetNumber]   ->Fill(  1, totalW);
			   hMassOutZevents[mc][jetNumber]->Fill(mll, totalW);
			 }
		       }
		     }
		   }		       
		 }
	       }
	     }

	     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	     //
	     // Main analisis
	     //
	     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	     
	     tree->Fill();
	     
	     hWTrigger[3]   ->Fill(1, totalW); 
	     hWeffTrigger[3]->Fill(1, efficiencyW);
	     
	     hNVtx[3]                           ->Fill(nvtx,      totalW); 
	     hPtLepton1TwoLeptonsLevel[3]       ->Fill(pt1,       totalW);
	     hPtLepton2TwoLeptonsLevel[3]       ->Fill(pt2,       totalW);
	     hPtDiLeptonTwoLeptonsLevel[3]      ->Fill(ptll,      totalW);
	     hMinvTwoLeptonsLevel[3]            ->Fill(mll,       totalW);
	     hMtTwoLeptonsLevel[3]              ->Fill(mth,       totalW);
	     hpfMetTwoLeptonsLevel[3]           ->Fill(pfType1Met,totalW);
	     htrkMetTwoLeptonsLevel[3]          ->Fill(trkMet,    totalW);
	     hpminMetTwoLeptonsLevel[3]         ->Fill(mpmet,     totalW);
	     hDeltaRLeptonsTwoLeptonsLevel[3]   ->Fill(drll,      totalW);
	     hDeltaPhiLeptonsTwoLeptonsLevel[3] ->Fill(dphill,    totalW);
	     hDPhiPtllJetTwoLeptonsLevel[3]     ->Fill(dphilljet, totalW);
	     hNjetsTwoLeptonsLevel[3]           ->Fill(njet,      totalW);
	     hNjetsPlot1TwoLeptonsLevel[3]      ->Fill(dRjet1,    totalW);
	     hNjetsPlot2TwoLeptonsLevel[3]      ->Fill(dRjet2,    totalW);
	     hSigMuNoHtTwoLeptonsLevel[3]       ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
	     hSigElNoHtTwoLeptonsLevel[3]       ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);
	     
	     if (nextra == 0) {
	       
	       hWExtraLepton[3]->Fill(1, totalW);
	       hWeffExtraLepton[3]->Fill(1, efficiencyW);
	       
	       if (mll > 12) {
		 
		 hWLowMinv[3]->Fill(1, totalW);
		 hWeffLowMinv[3]->Fill(1, efficiencyW);
		 
		 if (pfType1Met > 20 ) { // removed for differential xsec
		   
		   hWMetCut[3]->Fill(1, totalW);
		   hWeffMetCut[3]->Fill(1, efficiencyW);
		   
		   //zveto (in case of same flavour)
		   if ( (fabs(ZMASS - mll) > 15  && 
			 (metvar > 45 ) )        ||
			channel == 2             || 
			channel == 3             ){
		     
		     hWZVeto[3]->Fill(1, totalW);
		     hWeffZVeto[3]->Fill(1, efficiencyW);
		     
		     if (mpmet > 20){
		       
		       hWpMetCut[3]->Fill(1, totalW);
		       hWeffpMetCut[3]->Fill(1, efficiencyW);
		       
		       if (dphiv || channel == 2 || channel == 3) {
			 
			 hWDeltaPhiJet[3]->Fill(1, totalW);
			 hWeffDeltaPhiJet[3]->Fill(1, efficiencyW);
			 
			 if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) {
			   
			   hWPtll[3]->Fill(1, totalW);			    
			   hWeffPtll[3]->Fill(1, efficiencyW);			    
			   
			   hWnJets[3]->Fill(njet, totalW);
			   hWeffnJets[3]->Fill(njet, efficiencyW);
			   
			   hWnBtaggedJets[3]->Fill(nbjet, totalW);
			   hWeffnBtaggedJets[3]->Fill(nbjet, efficiencyW);
			   
			   hHt[3]->Fill(Ht,totalW);				    
			   /*
			     for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			     if (jetbin >= 3) jetbin = 2;
			     if(jetNumber == jetbin){
			     hHt[jetNumber]->Fill(Ht,totalW);				    
			     }
			     }
			   */
			   //b-veto
			   if (bveto_ip == 1 && nbjettche == 0)
			     hbVetoOld[3] -> Fill(1, totalW);              
			   
			   if (!bvetocsvv2ivfLoose)
			     hbVetoCsvv2ivfLoose[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfMedium)
			     hbVetoCsvv2ivfMedium[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfTight)
			     hbVetoCsvv2ivfTight[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfLoose && !bvetoMu)
			     hbVetoCsvv2ivfLooseAndMu[3] -> Fill(1, totalW);
			     //hHt[2]->Fill(Ht,totalW);				    
			     
			   if (!bvetocsvv2ivfRecommended && !bvetoMu)
			     hbVetoCsvv2ivfRecommendedAndMu[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfRecommended){
			     hbVetoCsvv2ivfRecommended[3] -> Fill(1, totalW);
			     
			     if (!bvetoMu){
			       hbVetoMu[3] -> Fill(1, totalW);
			     
			       hHt[2]->Fill(Ht,totalW);				    
			     

			       //bveto Ht 
			       if(Ht < 237){
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
				 hWnJetsBvetoAfterHt[3]     ->Fill(njet, efficiencyW);					
				 hSigMu[3]                 ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
				 hSigEl[3]                 ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);
		
				 //Control Region Selections - on Top of WW Cuts
				 if (drll > 1.5){
				   hpfMetCR[3]  -> Fill (pfType1Met, totalW);
				   htrkMetCR[3] -> Fill (trkMet,     totalW);
				   hmpMetCR[3]  -> Fill (mpmet,      totalW);
				 }
				 
				 //monoHiggs Selections
				 if (Mc < 100){
				   
				   hPtLepton1monoHLevel[0][3]      ->Fill(pt1,       totalW);
				   hPtLepton2monoHLevel[0][3]      ->Fill(pt2,       totalW);
				   hPtDiLeptonmonoHLevel[0][3]     ->Fill(ptll,      totalW);
				   hMinvmonoHLevel[0][3]           ->Fill(mll,       totalW);
				   hMtmonoHLevel[0][3]             ->Fill(mth,       totalW);
				   hMt1monoHLevel[0][3]            ->Fill(Mt1,       totalW);
				   hMt2monoHLevel[0][3]            ->Fill(Mt2,       totalW);
				   hpfMetmonoHLevel[0][3]          ->Fill(pfType1Met,totalW);
				   hpminMetmonoHLevel[0][3]        ->Fill(mpmet,     totalW);
				   hDeltaRLeptonsmonoHLevel[0][3]  ->Fill(drll,      totalW);
				   hDeltaPhiLeptonsmonoHLevel[0][3]->Fill(dphill,    totalW);
				   hDPhiPtllJetmonoHLevel[0][3]    ->Fill(dphilljet, totalW);
				   hDPhillMetmonoHLevel[0][3]      ->Fill(dphillmet, totalW);
				   hPtWWmonoHLevel[0][3]           ->Fill(ptWW,      totalW);
				   hMcmonoHLevel[0][3]             ->Fill(Mc,        totalW);		  
				   hTrkMetmonoHLevel[0][3]         ->Fill(trkMet,    totalW);
				   hHtmonoHLevel[0][3]             ->Fill(Ht,        totalW);			
				   hdphijetjetmonoHLevel[0][3]     ->Fill(dphijetjet,totalW);
				   hdetajetjetmonoHLevel[0][3]     ->Fill(detajetjet,totalW);
				   
				   if (drll < 1.5){
				     
				     hPtLepton1monoHLevel[1][3]      ->Fill(pt1,       totalW);
				     hPtLepton2monoHLevel[1][3]      ->Fill(pt2,       totalW);
				     hPtDiLeptonmonoHLevel[1][3]     ->Fill(ptll,      totalW);
				     hMinvmonoHLevel[1][3]           ->Fill(mll,       totalW);
				     hMtmonoHLevel[1][3]             ->Fill(mth,       totalW);
				     hMt1monoHLevel[1][3]            ->Fill(Mt1,       totalW);
				     hMt2monoHLevel[1][3]            ->Fill(Mt2,       totalW);
				     hpfMetmonoHLevel[1][3]          ->Fill(pfType1Met,totalW);
				     hpminMetmonoHLevel[1][3]        ->Fill(mpmet,     totalW);
				     hDeltaRLeptonsmonoHLevel[1][3]  ->Fill(drll,      totalW);
				     hDeltaPhiLeptonsmonoHLevel[1][3]->Fill(dphill,    totalW);
				     hDPhiPtllJetmonoHLevel[1][3]    ->Fill(dphilljet, totalW);
				     hDPhillMetmonoHLevel[1][3]      ->Fill(dphillmet, totalW);
				     hPtWWmonoHLevel[1][3]           ->Fill(ptWW,      totalW);
				     hMcmonoHLevel[1][3]             ->Fill(Mc,        totalW);		  
				     hTrkMetmonoHLevel[1][3]         ->Fill(trkMet,    totalW);
				     hHtmonoHLevel[1][3]             ->Fill(Ht,        totalW);			
				     hdphijetjetmonoHLevel[1][3]     ->Fill(dphijetjet,totalW);
				     hdetajetjetmonoHLevel[1][3]     ->Fill(detajetjet,totalW);
				     
				     if ( mpmet > 60 ){
				       
				       hPtLepton1monoHLevel[2][3]      ->Fill(pt1,       totalW);
				       hPtLepton2monoHLevel[2][3]      ->Fill(pt2,       totalW);
				       hPtDiLeptonmonoHLevel[2][3]     ->Fill(ptll,      totalW);
				       hMinvmonoHLevel[2][3]           ->Fill(mll,       totalW);
				       hMtmonoHLevel[2][3]             ->Fill(mth,       totalW);
				       hMt1monoHLevel[2][3]            ->Fill(Mt1,       totalW);
				       hMt2monoHLevel[2][3]            ->Fill(Mt2,       totalW);
				       hpfMetmonoHLevel[2][3]          ->Fill(pfType1Met,totalW);
				       hpminMetmonoHLevel[2][3]        ->Fill(mpmet,     totalW);
				       hDeltaRLeptonsmonoHLevel[2][3]  ->Fill(drll,      totalW);
				       hDeltaPhiLeptonsmonoHLevel[2][3]->Fill(dphill,    totalW);
				       hDPhiPtllJetmonoHLevel[2][3]    ->Fill(dphilljet, totalW);
				       hDPhillMetmonoHLevel[2][3]      ->Fill(dphillmet, totalW);
				       hPtWWmonoHLevel[2][3]           ->Fill(ptWW,      totalW);
				       hMcmonoHLevel[2][3]             ->Fill(Mc,        totalW);		  
				       hTrkMetmonoHLevel[2][3]         ->Fill(trkMet,    totalW);
				       hHtmonoHLevel[2][3]             ->Fill(Ht,        totalW);
				       hdphijetjetmonoHLevel[2][3]     ->Fill(dphijetjet,totalW);
				       hdetajetjetmonoHLevel[2][3]     ->Fill(detajetjet,totalW);
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
			   
	     for (int jetNumber = 0; jetNumber < 3 ; ++jetNumber){
	       if (jetbin >= 3) jetbin = 2;
	       if(jetNumber == jetbin){
		 
		 hWTrigger[jetNumber]   ->Fill(1, totalW); 
		 hWeffTrigger[jetNumber]->Fill(1, efficiencyW);
		 
		 hNVtx[jetNumber]                           ->Fill(nvtx,      totalW); 
		 hPtLepton1TwoLeptonsLevel[jetNumber]       ->Fill(pt1,       totalW);
		 hPtLepton2TwoLeptonsLevel[jetNumber]       ->Fill(pt2,       totalW);
		 hPtDiLeptonTwoLeptonsLevel[jetNumber]      ->Fill(ptll,      totalW);
		 hMinvTwoLeptonsLevel[jetNumber]            ->Fill(mll,       totalW);
		 hMtTwoLeptonsLevel[jetNumber]              ->Fill(mth,       totalW);
		 hpfMetTwoLeptonsLevel[jetNumber]           ->Fill(pfType1Met,totalW);
		 htrkMetTwoLeptonsLevel[jetNumber]          ->Fill(trkMet,    totalW);
		 hpminMetTwoLeptonsLevel[jetNumber]         ->Fill(mpmet,     totalW);
		 hDeltaRLeptonsTwoLeptonsLevel[jetNumber]   ->Fill(drll,      totalW);
		 hDeltaPhiLeptonsTwoLeptonsLevel[jetNumber] ->Fill(dphill,    totalW);
		 hDPhiPtllJetTwoLeptonsLevel[jetNumber]     ->Fill(dphilljet, totalW);
		 hNjetsTwoLeptonsLevel[jetNumber]           ->Fill(njet,      totalW);
		 hNjetsPlot1TwoLeptonsLevel[jetNumber]      ->Fill(dRjet1,    totalW);
		 hNjetsPlot2TwoLeptonsLevel[jetNumber]      ->Fill(dRjet2,    totalW);
		 hSigMuNoHtTwoLeptonsLevel[jetNumber]       ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
		 hSigElNoHtTwoLeptonsLevel[jetNumber]       ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);
		 
		 if (nextra == 0) {
		   
		   hWExtraLepton[jetNumber]->Fill(1, totalW);
		   hWeffExtraLepton[jetNumber]->Fill(1, efficiencyW);
		   
		   if (mll > 12) {
		     
		     hWLowMinv[jetNumber]->Fill(1, totalW);
		     hWeffLowMinv[jetNumber]->Fill(1, efficiencyW);
		     
		     if (pfType1Met > 20 ) { // removed for differential xsec
		       
		       hWMetCut[jetNumber]->Fill(1, totalW);
		       hWeffMetCut[jetNumber]->Fill(1, efficiencyW);
		       
		       //zveto (in case of same flavour)
		       if ( (fabs(ZMASS - mll) > 15  && 
			     (metvar > 45 ) )        ||
			    channel == 2             || 
			    channel == 3             ){
			 
			 hWZVeto[jetNumber]->Fill(1, totalW);
			 hWeffZVeto[jetNumber]->Fill(1, efficiencyW);
			 
			 if (mpmet > 20){
			   
			   hWpMetCut[jetNumber]->Fill(1, totalW);
			   hWeffpMetCut[jetNumber]->Fill(1, efficiencyW);
			   
			   if (dphiv || channel == 2 || channel == 3) {
			     
			     hWDeltaPhiJet[jetNumber]->Fill(1, totalW);
			     hWeffDeltaPhiJet[jetNumber]->Fill(1, efficiencyW);
			     
			     if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) {
			       
			       hWPtll[jetNumber]->Fill(1, totalW);			    
			       hWeffPtll[jetNumber]->Fill(1, efficiencyW);			    
			       
			       hWnJets[jetNumber]->Fill(njet, totalW);
			       hWeffnJets[jetNumber]->Fill(njet, efficiencyW);
			       
			       hWnBtaggedJets[jetNumber]->Fill(nbjet, totalW);
			       hWeffnBtaggedJets[jetNumber]->Fill(nbjet, efficiencyW);
			       
			       //b-veto
			       if (bveto_ip == 1 && nbjettche == 0)
				 hbVetoOld[jetNumber] -> Fill(1, totalW);              
			       
			       if (!bvetocsvv2ivfLoose && !bvetoMu)
				 hbVetoCsvv2ivfLooseAndMu[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfLoose)
				 hbVetoCsvv2ivfLoose[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfMedium)
				 hbVetoCsvv2ivfMedium[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfTight)
				 hbVetoCsvv2ivfTight[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfRecommended && !bvetoMu)
				 hbVetoCsvv2ivfRecommendedAndMu[jetNumber] -> Fill(1, totalW);
				 
			       if (!bvetocsvv2ivfRecommended){
				 hbVetoCsvv2ivfRecommended[jetNumber] -> Fill(1, totalW);
			       
				 if (!bvetoMu){
				   hbVetoMu[jetNumber] -> Fill(1, totalW);
				   
				   //bveto Ht  
				   if(Ht < 237){
				     
				     hPtLepton1WWLevel[jetNumber]      ->Fill(pt1,       totalW);
				     hPtLepton2WWLevel[jetNumber]      ->Fill(pt2,       totalW);
				     hPtDiLeptonWWLevel[jetNumber]     ->Fill(ptll,      totalW);
				     hMinvWWLevel[jetNumber]           ->Fill(mll,       totalW);
				     hMtWWLevel[jetNumber]             ->Fill(mth,       totalW);
				     hpfMetWWLevel[jetNumber]          ->Fill(pfType1Met,totalW);
				     htrkMetWWLevel[jetNumber]         ->Fill(trkMet,    totalW);
				     hpminMetWWLevel[jetNumber]        ->Fill(mpmet,     totalW);
				     hDeltaRLeptonsWWLevel[jetNumber]  ->Fill(drll,      totalW);
				     hDeltaPhiLeptonsWWLevel[jetNumber]->Fill(dphill,    totalW);
				     hDPhiPtllJetWWLevel[jetNumber]    ->Fill(dphilljet, totalW);
				     hSigMu[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
				     hSigEl[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);


				     //Control Region Selections - on Top of WW Cuts
				     if (drll > 1.5){
				       hpfMetCR[jetNumber]  -> Fill (pfType1Met, totalW);
				       htrkMetCR[jetNumber] -> Fill (trkMet,     totalW);
				       hmpMetCR[jetNumber]  -> Fill (mpmet,      totalW);
				     }
				     
				     //monoHiggs Selections
				     if (Mc < 100){
				       
				       hPtLepton1monoHLevel[0][jetNumber]      ->Fill(pt1,       totalW);
				       hPtLepton2monoHLevel[0][jetNumber]      ->Fill(pt2,       totalW);
				       hPtDiLeptonmonoHLevel[0][jetNumber]     ->Fill(ptll,      totalW);
				       hMinvmonoHLevel[0][jetNumber]           ->Fill(mll,       totalW);
				       hMtmonoHLevel[0][jetNumber]             ->Fill(mth,       totalW);
				       hMt1monoHLevel[0][jetNumber]            ->Fill(Mt1,       totalW);
				       hMt2monoHLevel[0][jetNumber]            ->Fill(Mt2,       totalW);
				       hpfMetmonoHLevel[0][jetNumber]          ->Fill(pfType1Met,totalW);
				       hpminMetmonoHLevel[0][jetNumber]        ->Fill(mpmet,     totalW);
				       hDeltaRLeptonsmonoHLevel[0][jetNumber]  ->Fill(drll,      totalW);
				       hDeltaPhiLeptonsmonoHLevel[0][jetNumber]->Fill(dphill,    totalW);
				       hDPhiPtllJetmonoHLevel[0][jetNumber]    ->Fill(dphilljet, totalW);
				       hDPhillMetmonoHLevel[0][jetNumber]      ->Fill(dphillmet, totalW);
				       hPtWWmonoHLevel[0][jetNumber]           ->Fill(ptWW,      totalW);
				       hMcmonoHLevel[0][jetNumber]             ->Fill(Mc,        totalW);		  
				       hTrkMetmonoHLevel[0][jetNumber]         ->Fill(trkMet,    totalW);
				       hHtmonoHLevel[0][jetNumber]             ->Fill(Ht,        totalW);			
				       hdphijetjetmonoHLevel[0][jetNumber]     ->Fill(dphijetjet,totalW);
				       hdetajetjetmonoHLevel[0][jetNumber]     ->Fill(detajetjet,totalW);
				       
				       if (drll < 1.5){
					 
					 hPtLepton1monoHLevel[1][jetNumber]      ->Fill(pt1,       totalW);
					 hPtLepton2monoHLevel[1][jetNumber]      ->Fill(pt2,       totalW);
					 hPtDiLeptonmonoHLevel[1][jetNumber]     ->Fill(ptll,      totalW);
					 hMinvmonoHLevel[1][jetNumber]           ->Fill(mll,       totalW);
					 hMtmonoHLevel[1][jetNumber]             ->Fill(mth,       totalW);
					 hMt1monoHLevel[1][jetNumber]            ->Fill(Mt1,       totalW);
					 hMt2monoHLevel[1][jetNumber]            ->Fill(Mt2,       totalW);
					 hpfMetmonoHLevel[1][jetNumber]          ->Fill(pfType1Met,totalW);
					 hpminMetmonoHLevel[1][jetNumber]        ->Fill(mpmet,     totalW);
					 hDeltaRLeptonsmonoHLevel[1][jetNumber]  ->Fill(drll,      totalW);
					 hDeltaPhiLeptonsmonoHLevel[1][jetNumber]->Fill(dphill,    totalW);
					 hDPhiPtllJetmonoHLevel[1][jetNumber]    ->Fill(dphilljet, totalW);
					 hDPhillMetmonoHLevel[1][jetNumber]      ->Fill(dphillmet, totalW);
					 hPtWWmonoHLevel[1][jetNumber]           ->Fill(ptWW,      totalW);
					 hMcmonoHLevel[1][jetNumber]             ->Fill(Mc,        totalW);		  
					 hTrkMetmonoHLevel[1][jetNumber]         ->Fill(trkMet,    totalW);
					 hHtmonoHLevel[1][jetNumber]             ->Fill(Ht,        totalW);			
					 hdphijetjetmonoHLevel[1][jetNumber]     ->Fill(dphijetjet,totalW);
					 hdetajetjetmonoHLevel[1][jetNumber]     ->Fill(detajetjet,totalW);
					 
					 if ( mpmet > 60 ){
					   
					   hPtLepton1monoHLevel[2][jetNumber]      ->Fill(pt1,       totalW);
					   hPtLepton2monoHLevel[2][jetNumber]      ->Fill(pt2,       totalW);
					   hPtDiLeptonmonoHLevel[2][jetNumber]     ->Fill(ptll,      totalW);
					   hMinvmonoHLevel[2][jetNumber]           ->Fill(mll,       totalW);
					   hMtmonoHLevel[2][jetNumber]             ->Fill(mth,       totalW);
					   hMt1monoHLevel[2][jetNumber]            ->Fill(Mt1,       totalW);
					   hMt2monoHLevel[2][jetNumber]            ->Fill(Mt2,       totalW);
					   hpfMetmonoHLevel[2][jetNumber]          ->Fill(pfType1Met,totalW);
					   hpminMetmonoHLevel[2][jetNumber]        ->Fill(mpmet,     totalW);
					   hDeltaRLeptonsmonoHLevel[2][jetNumber]  ->Fill(drll,      totalW);
					   hDeltaPhiLeptonsmonoHLevel[2][jetNumber]->Fill(dphill,    totalW);
					   hDPhiPtllJetmonoHLevel[2][jetNumber]    ->Fill(dphilljet, totalW);
					   hDPhillMetmonoHLevel[2][jetNumber]      ->Fill(dphillmet, totalW);
					   hPtWWmonoHLevel[2][jetNumber]           ->Fill(ptWW,      totalW);
					   hMcmonoHLevel[2][jetNumber]             ->Fill(Mc,        totalW);		  
					   hTrkMetmonoHLevel[2][jetNumber]         ->Fill(trkMet,    totalW);
					   hHtmonoHLevel[2][jetNumber]             ->Fill(Ht,        totalW);
					   hdphijetjetmonoHLevel[2][jetNumber]     ->Fill(dphijetjet,totalW);
					   hdetajetjetmonoHLevel[2][jetNumber]     ->Fill(detajetjet,totalW);
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
	     }
	   }

   // Define Normalization Factor for MC samples 
   //------------------------------------------------------------------------------
   
   // Define weights
   //------------------------------------------------------------------------------
   
   //float pileupweight = 1;
   
   //  if (!IsDATA)
   // pileupweight = fPUWeight->GetWeight(T_Event_nPU);
   
   //double factN = 1;
   
   //if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents;
   
   //factN = factN*pileupweight;
   //------------------------------------------------------------------------------
   
   // Init variables
   //------------------------------------------------------------------------------
   
} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//

   
   
void WWAnalysisSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master

  tree = FindOutput<TTree*>("nt");
  ///*** 1D histos ***/// 

  h_n_PV  = FindOutput<TH1F*>("h_N_PV");  
  
  //Z + Jets Data Driven histograms
  //----------------------------------------------------------------------------       
  char name[80];

  for (int j = 0; j < 4; ++j){
    for (size_t nC=0; nC<numberMetCuts; nC++) {
      sprintf(name,"hNinZevents%.1i%.1i",MetCut[nC],j);
      cout<<name<<endl;
      hNinZevents     [nC][j] = FindOutput<TH1F*>(name);
      sprintf(name,"hNoutZevents%.1i%.1i",MetCut[nC],j);
      hNoutZevents    [nC][j] = FindOutput<TH1F*>(name);
      sprintf(name,"hNinLooseZevents%.1i%.1i",MetCut[nC],j);
      hNinLooseZevents[nC][j] = FindOutput<TH1F*>(name);
      sprintf(name,"hMassInZevents%.1i%.1i",MetCut[nC],j);
      hMassInZevents  [nC][j] = FindOutput<TH1F*>(name);
      sprintf(name,"hMassOutZevents%.1i%.1i",MetCut[nC],j);
      hMassOutZevents [nC][j] = FindOutput<TH1F*>(name);
    }
  } 

  // Counting histograms
  //----------------------------------------------------------------------------

  for (Int_t qq = 0; qq < 4; ++qq){

    sprintf(name,"hWTrigger%.1i",qq);
    hWTrigger[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWMetCut%.1i",qq);
    hWMetCut[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hWLowMinv%.1i",qq);
    hWLowMinv[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWZVeto%.1i",qq);
    hWZVeto[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hWpMetCut%.1i",qq);
    hWpMetCut[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWJetVeto%.1i",qq);
    hWJetVeto[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWnJets%.1i",qq);
    hWnJets[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffnJets%.1i",qq);
    hWeffnJets[qq]    = FindOutput<TH1F*>(name);
    sprintf(name,"hWnBtaggedJets%.1i",qq);
    hWnBtaggedJets[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffnBtaggedJets%.1i",qq);
    hWeffnBtaggedJets[qq]  = FindOutput<TH1F*>(name);
    sprintf(name,"hWnJetsBveto%.1i",qq);
    hWnJetsBveto[qq]    = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffnJetsBveto%.1i",qq);
    hWeffnJetsBveto[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hNjetsTwoLeptonsLevel%.1i",qq);
    hNjetsTwoLeptonsLevel[qq]    = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffnJetsBvetoAfterHt%.1i",qq);
    hWeffnJetsBvetoAfterHt[qq] = FindOutput<TH1F*>(name);
    
    sprintf(name,"hWDeltaPhiJet%.1i",qq);
    hWDeltaPhiJet[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hWSoftMuVeto%.1i",qq);
    hWSoftMuVeto[qq]  = FindOutput<TH1F*>(name);
    sprintf(name,"hWExtraLepton%.1i",qq);
    hWExtraLepton[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hWPtll%.1i",qq);
    hWPtll[qq]        = FindOutput<TH1F*>(name);
    sprintf(name,"hWTopTagging%.1i",qq);
    hWTopTagging[qq]  = FindOutput<TH1F*>(name);
    
    sprintf(name,"hWeffTrigger%.1i",qq);
    hWeffTrigger[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffMetCut%.1i",qq);
    hWeffMetCut[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffLowMinv%.1i",qq);
    hWeffLowMinv[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffZVeto%.1i",qq);
    hWeffZVeto[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffpMetCut%.1i",qq);
    hWeffpMetCut[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffJetVeto%.1i",qq);
    hWeffJetVeto[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffDeltaPhiJet%.1i",qq);
    hWeffDeltaPhiJet[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffSoftMuVeto%.1i",qq);
    hWeffSoftMuVeto[qq]  = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffExtraLepton%.1i",qq);
    hWeffExtraLepton[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffPtll%.1i",qq);
    hWeffPtll[qq]        = FindOutput<TH1F*>(name);
    sprintf(name,"hWeffTopTagging%.1i",qq);
    hWeffTopTagging[qq]  = FindOutput<TH1F*>(name);
    
    sprintf(name,"hLooseIso%.1i",qq);
    hLooseIso[qq] = FindOutput<TH1F*>(name);
    
    // TwoLeptons level histograms
    //---------------------------------------------------------------------------
    sprintf(name,"hNVtx%.1i",qq);
    hNVtx[qq]                           = FindOutput<TH1F*>(name);
    sprintf(name,"hPtLepton1TwoLeptonsLevel%.1i",qq);
    hPtLepton1TwoLeptonsLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtLepton2TwoLeptonsLevel%.1i",qq);
    hPtLepton2TwoLeptonsLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hPtDiLeptonTwoLeptonsLevel%.1i",qq);
    hPtDiLeptonTwoLeptonsLevel[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hMinvTwoLeptonsLevel%.1i",qq);
    hMinvTwoLeptonsLevel[qq]            = FindOutput<TH1F*>(name);
    sprintf(name,"hMtTwoLeptonsLevel%.1i",qq);
    hMtTwoLeptonsLevel[qq]              = FindOutput<TH1F*>(name);
    sprintf(name,"hpfMetTwoLeptonsLevel%.1i",qq);
    hpfMetTwoLeptonsLevel[qq]           = FindOutput<TH1F*>(name);
    sprintf(name,"htrkMetTwoLeptonsLevel%.1i",qq);
    htrkMetTwoLeptonsLevel[qq]          = FindOutput<TH1F*>(name);
    sprintf(name,"hpminMetTwoLeptonsLevel%.1i",qq);
    hpminMetTwoLeptonsLevel[qq]         = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaRLeptonsTwoLeptonsLevel%.1i",qq);
    hDeltaRLeptonsTwoLeptonsLevel[qq]   = FindOutput<TH1F*>(name);
    sprintf(name,"hDeltaPhiLeptonsTwoLeptonsLevel%.1i",qq);
    hDeltaPhiLeptonsTwoLeptonsLevel[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hDPhiPtllJetTwoLeptonsLevel%.1i",qq);
    hDPhiPtllJetTwoLeptonsLevel[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hNjetsPlot1TwoLeptonsLevel%.1i",qq);
    hNjetsPlot1TwoLeptonsLevel[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hNjetsPlot2TwoLeptonsLevel%.1i",qq);
    hNjetsPlot2TwoLeptonsLevel[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hSigMuNoHtTwoLeptonsLevel%.1i",qq);
    hSigMuNoHtTwoLeptonsLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hSigElNoHtTwoLeptonsLevel%.1i",qq);
    hSigElNoHtTwoLeptonsLevel[qq]       = FindOutput<TH1F*>(name);
    sprintf(name,"hDxyTwoLeptonsLevel%.1i",qq);
    hDxyTwoLeptonsLevel[qq]             = FindOutput<TH1F*>(name);
    sprintf(name,"hDzTwoLeptonsLevel%.1i",qq);
    hDzTwoLeptonsLevel[qq]              = FindOutput<TH1F*>(name);
    
    sprintf(name,"hsoftMuPt%.1i",qq);
    hsoftMuPt[qq]                       = FindOutput<TH1F*>(name);
    sprintf(name,"hjetPt%.1i",qq);
    hjetPt[qq]                          = FindOutput<TH1F*>(name);
    
    // WW level histograms
    //----------------------------------------------------------------------------                                                    
    
    sprintf(name,"hWnJetsBvetoAfterHt%.1i",qq);
    hWnJetsBvetoAfterHt[qq] = FindOutput<TH1F*>(name);
    
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
    sprintf(name,"htrkMetWWLevel%.1i",qq);    
    htrkMetWWLevel[qq]          = FindOutput<TH1F*>(name);
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
    sprintf(name,"hHt%.1i",qq);
    hHt[qq]                     = FindOutput<TH1F*>(name);
    sprintf(name,"hHtAfter%.1i",qq);
    hHtAfter[qq]                = FindOutput<TH1F*>(name);
    
    // BVeto level histograms 
    //----------------------------------------------------------------------------
    
    sprintf(name,"hbVetoOld%.1i",qq);
    hbVetoOld[qq]                = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoMu%.1i",qq);
    hbVetoMu[qq]                 = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfLoose%.1i",qq);
    hbVetoCsvv2ivfLoose[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfMedium%.1i",qq);
    hbVetoCsvv2ivfMedium[qq]     = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfTight%.1i",qq);
    hbVetoCsvv2ivfTight[qq]      = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfLooseAndMu%.1i",qq);
    hbVetoCsvv2ivfLooseAndMu[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfRecommended%.1i",qq);
    hbVetoCsvv2ivfRecommended[qq] = FindOutput<TH1F*>(name);
    sprintf(name,"hbVetoCsvv2ivfRecommendedAndMu%.1i",qq);
    hbVetoCsvv2ivfRecommendedAndMu[qq] = FindOutput<TH1F*>(name);

    //Control Region Plots
    //----------------------------------------------------------------------------                
    
    sprintf(name,"hpfMetCR%.1i",qq);
    hpfMetCR[qq]  = FindOutput<TH1F*>(name);
    sprintf(name,"htrkMetCR%.1i",qq);
    htrkMetCR[qq] = FindOutput<TH1F*>(name); 
    sprintf(name,"hmpMetCR%.1i",qq);
    hmpMetCR[qq]  = FindOutput<TH1F*>(name); 

    // monoH level histograms     
    //----------------------------------------------------------------------------  
    
    for (Int_t jj=0; jj<4; jj++) {
      
      sprintf(name,"hPtLepton1monoHLevel%.1i%.1i",qq,jj);
      hPtLepton1monoHLevel[qq][jj]       = FindOutput<TH1F*>(name); 
      sprintf(name,"hPtLepton2monoHLevel%.1i%.1i",qq,jj);
      hPtLepton2monoHLevel[qq][jj]       = FindOutput<TH1F*>(name); 
      sprintf(name,"hPtDiLeptonmonoHLevel%.1i%.1i",qq,jj);
      hPtDiLeptonmonoHLevel[qq][jj]      = FindOutput<TH1F*>(name); 
      sprintf(name,"hMinvmonoHLevel%.1i%.1i",qq,jj);
      hMinvmonoHLevel[qq][jj]            = FindOutput<TH1F*>(name); 
      sprintf(name,"hMtmonoHLevel%.1i%.1i",qq,jj);
      hMtmonoHLevel[qq][jj]              = FindOutput<TH1F*>(name); 
      sprintf(name,"hMt1monoHLevel%.1i%.1i",qq,jj);
      hMt1monoHLevel[qq][jj]             = FindOutput<TH1F*>(name); 
      sprintf(name,"hMt2monoHLevel%.1i%.1i",qq,jj);
      hMt2monoHLevel[qq][jj]             = FindOutput<TH1F*>(name); 
      sprintf(name,"hpfMetmonoHLevel%.1i%.1i",qq,jj);
      hpfMetmonoHLevel[qq][jj]           = FindOutput<TH1F*>(name); 
      sprintf(name,"hpminMetmonoHLevel%.1i%.1i",qq,jj);
      hpminMetmonoHLevel[qq][jj]         = FindOutput<TH1F*>(name); 
      sprintf(name,"hDeltaRLeptonsmonoHLevel%.1i%.1i",qq,jj);
      hDeltaRLeptonsmonoHLevel[qq][jj]   = FindOutput<TH1F*>(name); 
      sprintf(name,"hDeltaPhiLeptonsmonoHLevel%.1i%.1i",qq,jj);
      hDeltaPhiLeptonsmonoHLevel[qq][jj] = FindOutput<TH1F*>(name); 
      sprintf(name,"hDPhiPtllJetmonoHLevel%.1i%.1i",qq,jj);
      hDPhiPtllJetmonoHLevel[qq][jj]     = FindOutput<TH1F*>(name); 
      sprintf(name,"hDPhillMetmonoHLevel%.1i%.1i",qq,jj);
      hDPhillMetmonoHLevel[qq][jj]       = FindOutput<TH1F*>(name); 
      sprintf(name,"hPtWWmonoHLevel%.1i%.1i",qq,jj);
      hPtWWmonoHLevel[qq][jj]            = FindOutput<TH1F*>(name); 
      sprintf(name,"hMcmonoHLevel%.1i%.1i",qq,jj);
      hMcmonoHLevel[qq][jj]              = FindOutput<TH1F*>(name); 
      sprintf(name,"hTrkMetmonoHLevel%.1i%.1i",qq,jj);
      hTrkMetmonoHLevel[qq][jj]          = FindOutput<TH1F*>(name); 
      sprintf(name,"hHtmonoHLevel%.1i%.1i",qq,jj);
      hHtmonoHLevel[qq][jj]              = FindOutput<TH1F*>(name); 
      sprintf(name,"hdphijetjetmonoHLevel%.1i%.1i",qq,jj);
      hdphijetjetmonoHLevel[qq][jj]      = FindOutput<TH1F*>(name); 
      sprintf(name,"hdetajetjetmonoHLevel%.1i%.1i",qq,jj);
      hdetajetjetmonoHLevel[qq][jj]      = FindOutput<TH1F*>(name); 
    }
  }
}

//------------------------------------------------------------------------------
// IsTightLepton
//
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
// egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium
//------------------------------------------------------------------------------
bool WWAnalysisSelector::IsTightLepton(int k, TString _MuonID_)
{
  bool is_tight_lepton = false;

  // Muon ID
  if (fabs(std_vector_lepton_flavour.at(k)) == 13){
    
    if (_MuonID_ == "MediumID"){
      if (std_vector_lepton_isMediumMuon.at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "MediumIDTighterIP"){
      if (std_vector_lepton_isMediumMuon.at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy.at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz.at(k)) < 0.1)
	    is_tight_lepton = true;
    }
    
    else if (_MuonID_ == "TightID" ){
      if (std_vector_lepton_isTightMuon.at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "TightIDTighterIP"){
      if (std_vector_lepton_isTightMuon.at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy.at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz.at(k)) < 0.1)
	    is_tight_lepton = true;
    }
  }

  else if (fabs(std_vector_lepton_flavour.at(k)) == 11){
    is_tight_lepton = std_vector_lepton_eleIdTight.at(k);
  }

  /*  
  // Electron cut based medium ID plus tighter trigger cuts
  else if (fabs(std_vector_lepton_flavour.at(k)) == 11){
    //ECAL Barrel
    if( fabs(std_vector_lepton_eta.at(k)  <= 1.479) &&
        std_vector_electron_hOverE.at(k)  < 0.08    && 
	std_vector_electron_dEtaIn.at(k)  < 0.01    &&
	std_vector_electron_dPhiIn.at(k)  < 0.04    &&
	std_vector_electron_ooEmooP.at(k) < 0.01   ){
      
      is_tight_lepton = std_vector_lepton_eleIdMedium.at(k);
    }
  //ECAL Endcap
    if( fabs(std_vector_lepton_eta.at(k) > 1.479) && 
        std_vector_electron_hOverE.at(k) < 0.08   && 
	std_vector_electron_dEtaIn.at(k) < 0.01   &&
	std_vector_electron_dPhiIn.at(k) < 0.08   &&
	std_vector_electron_ooEmooP.at(k) < 0.01  ){
      
      is_tight_lepton = std_vector_lepton_eleIdMedium.at(k);
    }
  }
  */
  /*   
  float aeta = fabs(std_vector_electron_scEta.at(k));
      
      if (aeta <= 1.479)
	{
	  if (fabs(std_vector_electron_deltaEtaIn.at(k)) < 0.008925 &&
	      fabs(std_vector_electron_deltaPhiIn.at(k)) < 0.035973 &&
	      std_vector_electron_sigmaIetaIeta.at(k)    < 0.009996 &&
	      std_vector_electron_HoE.at(k)              < 0.050537 &&
	      fabs(std_vector_electron_d0.at(k))         < 0.012235 &&
	      fabs(std_vector_electron_dz.at(k))         < 0.042020 &&
	      fabs(std_vector_electron_ooEooP.at(k))     < 0.091942 &&
	      ElectronIsolation(k)                        < 0.107587 &&
	      !std_vector_electron_passConversion.at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
      else if (aeta > 1.479 && aeta < 2.5)
	{
	  if (fabs(std_vector_electron_deltaEtaIn.at(k)) < 0.007429 &&
	      fabs(std_vector_electron_deltaPhiIn.at(k)) < 0.067879 &&
	      std_vector_electron_sigmaIetaIeta.at(k)    < 0.030135 &&
	      std_vector_electron_HoE.at(k)              < 0.086782 &&
	      fabs(std_vector_electron_d0.at(k))         < 0.036719 &&
	      fabs(std_vector_electron_dz.at(k))         < 0.138142 &&
	      fabs(std_vector_electron_ooEooP.at(k))     < 0.100683 &&
	      ElectronIsolation(k)                        < 0.113254 &&
	      !std_vector_electron_passConversion.at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
    }
  */
  return is_tight_lepton;

}


//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float WWAnalysisSelector::MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt.at(k);
  float id = std_vector_lepton_flavour.at(k);

  float relative_isolation = -999;

  if (fabs(id) != 13) return relative_isolation;
  
  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso.at(k) +
	      std_vector_lepton_neutralHadronIso.at(k) -
	      0.5*std_vector_lepton_sumPUPt.at(k)));
  
  /*
  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) + 
    std_vector_lepton_neutralHadronIso.at(k) + 
    std_vector_lepton_photonIso.at(k);
  */
  relative_isolation =  relative_isolation / pt;

  return relative_isolation;
}


//------------------------------------------------------------------------------
// ElectronIsolation
//------------------------------------------------------------------------------
float WWAnalysisSelector::ElectronIsolation(int k)
{
  float pt = std_vector_lepton_pt.at(k);
  float id = std_vector_lepton_flavour.at(k);

  float relative_isolation = -999;

  if (fabs(id) != 11) return relative_isolation;

  
  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso.at(k) +
	      std_vector_lepton_neutralHadronIso.at(k) -
	      jetRho*std_vector_electron_effectiveArea.at(k)));
  
  /*
  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) + 
    std_vector_lepton_neutralHadronIso.at(k) + 
    std_vector_lepton_photonIso.at(k);
  */
  relative_isolation = relative_isolation / pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// IsIsolatedLepton
//------------------------------------------------------------------------------
bool WWAnalysisSelector::IsIsolatedLepton(int k)
{
  float id = std_vector_lepton_flavour.at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}
