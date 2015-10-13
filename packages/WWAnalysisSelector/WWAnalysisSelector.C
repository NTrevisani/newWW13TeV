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
  
  hWTrigger     = CreateH1F("hWTrigger",                        "", 10, 0, 10);
  hWMetCut      = CreateH1F("hWMetCut",                         "", 10, 0, 10);
  hWLowMinv     = CreateH1F("hWLowMinv",                        "", 10, 0, 10);
  hWZVeto       = CreateH1F("hWZVeto",                          "", 10, 0, 10);
  hWpMetCut     = CreateH1F("hWpMetCut",                        "", 10, 0, 10);
  hWJetVeto     = CreateH1F("hWJetVeto",                        "", 10, 0, 10);
  hWnJets       = CreateH1F("hWnJets",                          "", 10, 0, 10);
  hWeffnJets    = CreateH1F("hWeffnJets",                       "", 10, 0, 10);
  hWnBtaggedJets     = CreateH1F("hWnBtaggedJets",              "", 10, 0, 10);
  hWeffnBtaggedJets  = CreateH1F("hWeffnBtaggedJets",           "", 10, 0, 10);
  hWnJetsBveto    = CreateH1F("hWnJetsBveto",                   "", 10, 0, 10);
  hWeffnJetsBveto = CreateH1F("hWeffnJetsBveto",                "", 10, 0, 10);
  hNjetsTwoLeptonsLevel    = CreateH1F("hNjetsTwoLeptonsLevel", "", 10, 0, 10);
  hWeffnJetsBvetoAfterHt = CreateH1F("hWeffnJetsBvetoAfterHt",  "", 10, 0, 10);
  
  hWDeltaPhiJet = CreateH1F("hWDeltaPhiJet", "", 10, 0, 10);
  hWSoftMuVeto  = CreateH1F("hWSoftMuVeto",  "", 10, 0, 10);
  hWExtraLepton = CreateH1F("hWExtraLepton", "", 10, 0, 10);
  hWPtll        = CreateH1F("hWPtll",        "", 10, 0, 10);
  hWTopTagging  = CreateH1F("hWTopTagging",  "", 10, 0, 10);

  hWeffTrigger     = CreateH1F("hWeffTrigger",     "", 10, 0, 10);
  hWeffMetCut      = CreateH1F("hWeffMetCut",      "", 10, 0, 10);
  hWeffLowMinv     = CreateH1F("hWeffLowMinv",     "", 10, 0, 10);
  hWeffZVeto       = CreateH1F("hWeffZVeto",       "", 10, 0, 10);
  hWeffpMetCut     = CreateH1F("hWeffpMetCut",     "", 10, 0, 10);
  hWeffJetVeto     = CreateH1F("hWeffJetVeto",     "", 10, 0, 10);
  hWeffDeltaPhiJet = CreateH1F("hWeffDeltaPhiJet", "", 10, 0, 10);
  hWeffSoftMuVeto  = CreateH1F("hWeffSoftMuVeto",  "", 10, 0, 10);
  hWeffExtraLepton = CreateH1F("hWeffExtraLepton", "", 10, 0, 10);
  hWeffPtll        = CreateH1F("hWeffPtll",        "", 10, 0, 10);
  hWeffTopTagging  = CreateH1F("hWeffTopTagging",  "", 10, 0, 10);
  
  hLooseIso = CreateH1F("hLooseIso",  "", 100, 0, 10);

  // TwoLeptons level histograms   
  //----------------------------------------------------------------------------

  hNVtx                           = CreateH1F("hNVtx",                           "",  100,0., 100);  
  hPtLepton1TwoLeptonsLevel       = CreateH1F("hPtLepton1TwoLeptonsLevel",       "", 3000,0.,3000);
  hPtLepton2TwoLeptonsLevel       = CreateH1F("hPtLepton2TwoLeptonsLevel",       "", 3000,0.,3000);
  hPtDiLeptonTwoLeptonsLevel      = CreateH1F("hPtDiLeptonTwoLeptonsLevel",      "", 3000,0.,3000);
  hMinvTwoLeptonsLevel            = CreateH1F("hMinvTwoLeptonsLevel",            "", 3000,0.,3000);
  hMtTwoLeptonsLevel              = CreateH1F("hMtTwoLeptonsLevel",              "", 3000,0.,3000);
  hpfMetTwoLeptonsLevel           = CreateH1F("hpfMetTwoLeptonsLevel",           "", 3000,0.,3000);
  htrkMetTwoLeptonsLevel          = CreateH1F("htrkMetTwoLeptonsLevel",          "", 3000,0.,3000);
  hpminMetTwoLeptonsLevel         = CreateH1F("hpminMetTwoLeptonsLevel",         "", 3000,0.,3000);
  hDeltaRLeptonsTwoLeptonsLevel   = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",  50, 0.,   5);
  hDeltaPhiLeptonsTwoLeptonsLevel = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",  32, 0., 3.2);
  hDPhiPtllJetTwoLeptonsLevel     = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",     "",  32, 0., 3.2);
  hNjetsPlot1TwoLeptonsLevel      = CreateH1F("hNjetsPlot1TwoLeptonsLevel",      "",  10, 0.,  10);
  hNjetsPlot2TwoLeptonsLevel      = CreateH1F("hNjetsPlot2TwoLeptonsLevel",      "",  10, 0.,  10);
  hSigMuNoHtTwoLeptonsLevel       = CreateH1F("hSigMuNoHtTwoLeptonsLevel",       "",  10, 0.,  10);
  hSigElNoHtTwoLeptonsLevel       = CreateH1F("hSigElNoHtTwoLeptonsLevel",       "",  10, 0.,  10);
  hDxyTwoLeptonsLevel             = CreateH1F("hDxyTwoLeptonsLevel",             "", 1000,-0.1,0.1);
  hDzTwoLeptonsLevel              = CreateH1F("hDzTwoLeptonsLevel",              "", 1000,-0.5,0.5);

  hsoftMuPt                       = CreateH1F("hsoftMuPt",                       "", 3000,0.,3000);
  hjetPt                          = CreateH1F("hjetPt",                          "", 3000,0.,3000);

  // WW level histograms     
  //----------------------------------------------------------------------------  
  
  hWnJetsBvetoAfterHt = CreateH1F("hWnJetsBvetoAfterHt",  "", 10, 0, 10); 

  for (Int_t nC=0; nC<4; nC++) {

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
  
  //NNLL corrections
  if (_Signal == "WW50")
    baseW = baseW * nllW;

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
  
  //building Ht
  Ht = std_vector_lepton_pt.at(0) + std_vector_lepton_pt.at(1) + pfType1Met;
  
  if(njet > 10) njet = 10;
  for (int i = 0; i < njet; ++i)
    if(std_vector_jet_pt.at(i) > 0)
      Ht += std_vector_jet_pt.at(i);

  //defining nextra
  nextra = 0;
  for (unsigned int i = 0; i < std_vector_lepton_pt.size(); ++i)
    if(std_vector_lepton_pt.at(i) > 10)
      ++nextra;
  nextra = nextra - 2;
       
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

   for (unsigned int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.423)
       bvetocsvv2ivfLoose = true;

   //Building b-veto csvv2ivf Medium (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfMedium = false;

   for (unsigned int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.814)
       bvetocsvv2ivfMedium = true;

   //Building b-veto csvv2ivf Tight (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfTight = false;

   for (unsigned int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.941)
       bvetocsvv2ivfTight = true;

   //Building b-veto csvv2ivf Recommended (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfRecommended = false;

   for (unsigned int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.605)
       bvetocsvv2ivfRecommended = true;

   //Building a Soft Muon B-Veto (Quite Rudely)
   bool bvetoMu = false;
   
   for (unsigned int i = 0; i < std_vector_jet_softMuPt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 10 && std_vector_jet_pt.at(i) < 30)
       if (std_vector_jet_softMuPt.at(i) > 3){
	 hsoftMuPt->Fill(std_vector_jet_softMuPt.at(i), totalW);
	 hjetPt   ->Fill(std_vector_jet_pt.at(i),       totalW);
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

   // The selection begins here
   //--------------------------------------------------------------------------
   if (trigger == 1)
     //if (std_vector_lepton_pt.at(0) > 20)
     //if (std_vector_lepton_pt.at(1) > 20) 
     if (a != 100)
       if (b != 100)
	 if ((_SameSign == "SS" && ch1*ch2 > 0) || (_SameSign == "OS" && ch1*ch2 < 0))
	   if ( (SelectedChannel == -1)                                     || 
		(channel == SelectedChannel)                                || 
		(SelectedChannel == 4 && (channel == 2 || channel == 3) )   || 
		(SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
		){
	     /*	    
		    if (IsTightLepton(0,_MuonID) && !IsTightLepton(1,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(1), totalW);
		    if (IsTightLepton(1,_MuonID) && !IsTightLepton(0,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(0), totalW);
	     */
	     
	     //lepton a
	     hDxyTwoLeptonsLevel -> Fill(std_vector_lepton_BestTrackdxy.at(a),totalW);
	     hDzTwoLeptonsLevel  -> Fill(std_vector_lepton_BestTrackdz.at(a) ,totalW);
	     
	     //lepton b
	     hDxyTwoLeptonsLevel -> Fill(std_vector_lepton_BestTrackdxy.at(b),totalW);
	     hDzTwoLeptonsLevel  -> Fill(std_vector_lepton_BestTrackdz.at(b) ,totalW);
	     
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
	     
	     hWTrigger   ->Fill(1, totalW); 
	     hWeffTrigger->Fill(1, efficiencyW);
	     
	     hNVtx                          ->Fill(nvtx,      totalW); 
	     hPtLepton1TwoLeptonsLevel      ->Fill(pt1,       totalW);
	     hPtLepton2TwoLeptonsLevel      ->Fill(pt2,       totalW);
	     hPtDiLeptonTwoLeptonsLevel     ->Fill(ptll,      totalW);
	     hMinvTwoLeptonsLevel           ->Fill(mll,       totalW);
	     hMtTwoLeptonsLevel             ->Fill(mth,       totalW);
	     hpfMetTwoLeptonsLevel          ->Fill(pfType1Met,totalW);
	     htrkMetTwoLeptonsLevel         ->Fill(trkMet,    totalW);
	     hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
	     hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
	     hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
	     hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
	     hNjetsTwoLeptonsLevel          ->Fill(njet,      totalW);
	     hNjetsPlot1TwoLeptonsLevel     ->Fill(dRjet1,    totalW);
	     hNjetsPlot2TwoLeptonsLevel     ->Fill(dRjet2,    totalW);
	     hSigMuNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
	     hSigElNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);
	     
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
		   if ( (fabs(ZMASS - mll) > 15  && 
			 (metvar > 45 ) )        ||
			channel == 2             || 
			channel == 3             ){
		     
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
			   
			   if (!bvetoMu)  
			     hbVetoMu[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfLoose)
			     hbVetoCsvv2ivfLoose[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfMedium)
			     hbVetoCsvv2ivfMedium[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfTight)
			     hbVetoCsvv2ivfTight[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfLoose && !bvetoMu)
			     hbVetoCsvv2ivfLooseAndMu[3] -> Fill(1, totalW);
			     //hHt[2]->Fill(Ht,totalW);				    
			     
			   if (!bvetocsvv2ivfRecommended)
			     hbVetoCsvv2ivfRecommended[3] -> Fill(1, totalW);
			   
			   if (!bvetocsvv2ivfRecommended && !bvetoMu){
			     hbVetoCsvv2ivfRecommendedAndMu[3] -> Fill(1, totalW);
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
			       hWnJetsBvetoAfterHt       ->Fill(njet, efficiencyW);					
			       hSigMu[3]                 ->Fill(std_vector_lepton_muSIP3D.at(a),totalW);
			       hSigEl[3]                 ->Fill(std_vector_lepton_elSIP3D.at(a),totalW);
			     }
			   }				     
			   
			   for (int jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			     if (jetbin >= 3) jetbin = 2;
			     if(jetNumber == jetbin){
			       
			       //b-veto
			       if (bveto_ip == 1 && nbjettche == 0)
				 hbVetoOld[jetNumber] -> Fill(1, totalW);              
			       
			       if (!bvetoMu)
				 hbVetoMu[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfLoose && !bvetoMu)
				 hbVetoCsvv2ivfLooseAndMu[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfLoose)
				 hbVetoCsvv2ivfLoose[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfMedium)
				 hbVetoCsvv2ivfMedium[jetNumber] -> Fill(1, totalW);
			       
			       if (!bvetocsvv2ivfTight)
				 hbVetoCsvv2ivfTight[jetNumber] -> Fill(1, totalW);
				 
				 if (!bvetocsvv2ivfRecommended)
				   hbVetoCsvv2ivfRecommended[jetNumber] -> Fill(1, totalW);
				 
				 if (!bvetocsvv2ivfRecommended && !bvetoMu){
				   hbVetoCsvv2ivfRecommendedAndMu[jetNumber] -> Fill(1, totalW);
			     
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
   
   double factN = 1;
   
   if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents;
   
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

  hWTrigger     = FindOutput<TH1F*>("hWTrigger");                        
  hWMetCut      = FindOutput<TH1F*>("hWMetCut");                         
  hWLowMinv     = FindOutput<TH1F*>("hWLowMinv");                        
  hWZVeto       = FindOutput<TH1F*>("hWZVeto");                          
  hWpMetCut     = FindOutput<TH1F*>("hWpMetCut");                        
  hWJetVeto     = FindOutput<TH1F*>("hWJetVeto");                        
  hWnJets       = FindOutput<TH1F*>("hWnJets");                          
  hWeffnJets    = FindOutput<TH1F*>("hWeffnJets");                       
  hWnBtaggedJets     = FindOutput<TH1F*>("hWnBtaggedJets");              
  hWeffnBtaggedJets  = FindOutput<TH1F*>("hWeffnBtaggedJets");           
  hWnJetsBveto    = FindOutput<TH1F*>("hWnJetsBveto");                   
  hWeffnJetsBveto = FindOutput<TH1F*>("hWeffnJetsBveto");                
  hNjetsTwoLeptonsLevel    = FindOutput<TH1F*>("hNjetsTwoLeptonsLevel"); 
  hWeffnJetsBvetoAfterHt = FindOutput<TH1F*>("hWeffnJetsBvetoAfterHt");  
  
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
  
  // TwoLeptons level histograms
  //---------------------------------------------------------------------------
  hNVtx                           = FindOutput<TH1F*>("hNVtx");   
  hPtLepton1TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton1TwoLeptonsLevel");       
  hPtLepton2TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton2TwoLeptonsLevel");       
  hPtDiLeptonTwoLeptonsLevel      = FindOutput<TH1F*>("hPtDiLeptonTwoLeptonsLevel");      
  hMinvTwoLeptonsLevel            = FindOutput<TH1F*>("hMinvTwoLeptonsLevel");            
  hMtTwoLeptonsLevel              = FindOutput<TH1F*>("hMtTwoLeptonsLevel");              
  hpfMetTwoLeptonsLevel           = FindOutput<TH1F*>("hpfMetTwoLeptonsLevel");           
  htrkMetTwoLeptonsLevel          = FindOutput<TH1F*>("htrkMetTwoLeptonsLevel");          
  hpminMetTwoLeptonsLevel         = FindOutput<TH1F*>("hpminMetTwoLeptonsLevel");         
  hDeltaRLeptonsTwoLeptonsLevel   = FindOutput<TH1F*>("hDeltaRLeptonsTwoLeptonsLevel");   
  hDeltaPhiLeptonsTwoLeptonsLevel = FindOutput<TH1F*>("hDeltaPhiLeptonsTwoLeptonsLevel"); 
  hDPhiPtllJetTwoLeptonsLevel     = FindOutput<TH1F*>("hDPhiPtllJetTwoLeptonsLevel");     
  hNjetsPlot1TwoLeptonsLevel      = FindOutput<TH1F*>("hNjetsPlot1TwoLeptonsLevel");      
  hNjetsPlot2TwoLeptonsLevel      = FindOutput<TH1F*>("hNjetsPlot2TwoLeptonsLevel");      
  hSigMuNoHtTwoLeptonsLevel       = FindOutput<TH1F*>("hSigMuNoHtTwoLeptonsLevel");       
  hSigElNoHtTwoLeptonsLevel       = FindOutput<TH1F*>("hSigElNoHtTwoLeptonsLevel");       
  hDxyTwoLeptonsLevel             = FindOutput<TH1F*>("hDxyTwoLeptonsLevel");             
  hDzTwoLeptonsLevel              = FindOutput<TH1F*>("hDzTwoLeptonsLevel");              

  hsoftMuPt                       = FindOutput<TH1F*>("hsoftMuPt");                       
  hjetPt                          = FindOutput<TH1F*>("hjetPt");                          

  // WW level histograms
  //----------------------------------------------------------------------------                                                    

  hWnJetsBvetoAfterHt = FindOutput<TH1F*>("hWnJetsBvetoAfterHt");

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
  
  // Electron cut based medium ID plus tighter trigger cuts
  else if (fabs(std_vector_lepton_flavour.at(k)) == 11){
    //ECAL Barrel
    if( fabs(std_vector_lepton_eta.at(k) < 0.8)   &&
        std_vector_electron_hOverE.at(k) < 0.08   && 
	std_vector_electron_dEtaIn.at(k) < 0.01   &&
	std_vector_electron_dPhiIn.at(k) < 0.04   &&
	std_vector_electron_ooEmooP.at(k) < 0.01  ){
      
      is_tight_lepton = std_vector_lepton_eleIdMedium.at(k);
    }
  //ECAL Endcap
    if( fabs(std_vector_lepton_eta.at(k) > 1.2 && fabs(std_vector_lepton_eta.at(k)) < 2.4) &&
        std_vector_electron_hOverE.at(k) < 0.08   && 
	std_vector_electron_dEtaIn.at(k) < 0.01   &&
	std_vector_electron_dPhiIn.at(k) < 0.08   &&
	std_vector_electron_ooEmooP.at(k) < 0.01  ){
      
      is_tight_lepton = std_vector_lepton_eleIdMedium.at(k);
    }
  }
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
