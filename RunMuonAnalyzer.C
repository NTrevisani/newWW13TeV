//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                     WW ANALYSIS WITH PAF                                    /////////////
/////////////                                                                                             /////////////
/////////////                                    Nicolò Trevisani (IFCA)                                  /////////////
/////////////                                          Jun 2016                                           /////////////
/////////////                                                                                             /////////////
/////////////                              -> Adjust to a 120 width window <-                             /////////////
/////////////                                                                                             /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TSelector.h"

void RunMuonAnalyzer(TString data,
		     TString flavorChannel,
		     TString sameSign,
		     TString PAFMode,
		     Float_t luminosity,  //In fb-1
		     TString muonID
		     ) {
 
  //  gSystem->Load("libPAF.so");   
  TString path25      = "/gpfs/csic_projects/tier3data/LatinosSkims/MC_Spring15/25ns_August_PU/";
  TString path50      = "/gpfs/csic_projects/tier3data/LatinosSkims/MC_Spring15/50ns_August_PU/";
  TString pathOld     = "/gpfs/csic_projects/cms/trevisanin/newLatino/";
  TString pathData    = "/gpfs/csic_projects/tier3data/LatinosSkims/Data13TeV/";

  TString outPath     = "/gpfs/csic_users/trevisanin/newWW13TeV/rootFiles/" + flavorChannel + "/" + muonID + "/";

  gSystem -> Exec("mkdir /gpfs/csic_users/trevisanin/newWW13TeV/rootFiles");
  gSystem -> Exec("mkdir /gpfs/csic_users/trevisanin/newWW13TeV/rootFiles/" + flavorChannel);
  gSystem -> Exec("mkdir /gpfs/csic_users/trevisanin/newWW13TeV/rootFiles/" + flavorChannel + "/" + muonID);
  gSystem -> Exec("mkdir /gpfs/csic_users/trevisanin/newWW13TeV/rootFiles/" + flavorChannel + "/" + muonID + "/25ns/");
  gSystem -> Exec("mkdir /gpfs/csic_users/trevisanin/newWW13TeV/rootFiles/" + flavorChannel + "/" + muonID + "/50ns/");

  // Manual Input Parameters
  bool     debug            = true;    //For verbose while debugging
  int      nEventsToProcess = -1;      //Number of events to be processed (-1 = all)
  bool     doReport         = false;   //Count events and print final report

  // Automatic Input Parameters (don't touch)
  bool     isdata             = false;
  int      whichRun           = 1;
  int      nEventsInTheSample = 1; //before skimming
  float    xSection           = 1.;
  
  //---------------------------------
  // INITIALISE PAF PROJECT
  //---------------------------------

  // Create Project in Sequential Environment mode
  if ( PAFMode == "Sequential" || PAFMode == "sequential")
    PAFProject* myProject = new PAFProject( new PAFSequentialEnvironment() );

  // Create Project in Cluster Environment mode
  else if ( PAFMode == "Cluster" || PAFMode == "cluster" )
  PAFProject* myProject = new PAFProject( new PAFPROOFClusterEnvironment(10,20) );

  // Create Project in Lite Environment mode
  else if ( PAFMode == "Lite" || PAFMode == "lite" )
  PAFProject* myProject = new PAFProject( new PAFPROOFLiteEnvironment(8) );

  // Create Project in PoD Environment mode
  else if ( PAFMode == "PoD" || PAFMode == "pod" || PAFMode == "POD" )
  PAFProject* myProject = new PAFProject( new PAFPoDEnvironment() );

  else {
    cout<<"Please select a valid PAF operating mode: Sequential or Cluster"<<endl;
    return; 
  }

  // Set default name and subdirectory of Trees
  myProject->SetDefaultTreeName("latino");

  ////////////////////////////////////////////////
  // ADD SAMPLES

  //25ns
  //*******************************************************************************

  if (data=="WW25") {

    myProject->AddDataFile(path25 + "latino_WWTo2L2Nu.root");

    isdata             = false;
    nEventsInTheSample = 538089; 
    xSection           = 12.461;
    whichRun           = 2; 
  }

  else if (data=="WJets25") {

    myProject->AddDataFile(path25 + "latino_WJetsToLNu.root");
    
    isdata             = false;
    nEventsInTheSample = 10370741; 
    xSection           = 60781.5;
    whichRun           = 2; 
  }

  else if (data=="ZZ25") {

    myProject->AddDataFile(path25 + "latino_ZZ.root");
    
    isdata             = false;
    nEventsInTheSample = 206249; 
    xSection           = 31.8;
    whichRun           = 2; 
  }

  else if (data=="singleTop25") {

    myProject->AddDataFile(path25 + "latino_topTchannelAntitop.root");

    isdata             = false;
    nEventsInTheSample = 1112252;
    xSection           = ;
    whichRun           = 2;
  }

  else if (data=="HWW25") {

    myProject->AddDataFile(path25 + "latino_GluGluHToWWTo2L2Nu_M125.root");

    isdata             = false;
    nEventsInTheSample = 475039 + 473153;
    xSection           = 71.6;
    whichRun           = 2;
  }

  else if (data=="TW25") {

    myProject->AddDataFile(path25 + "latino_ST_tW_antitop.root");
    myProject->AddDataFile(path25 + "latino_ST_tW_top.root");

    isdata             = false;
    nEventsInTheSample = 372091;
    xSection           = 0.9913;
    whichRun           = 2;
  }

  else if (data=="DY25") {

    myProject->AddDataFile(path25 + "latino_DYJetsToLL_M-10to50.root");

    isdata             = false;
    nEventsInTheSample = 4081613;
    xSection           = 18610;
    whichRun           = 2;
  }

  //50ns
  //*******************************************************************************

  else if (data=="Data201550") {
    /*
    myProject->AddDataFile(pathData + "50ns/latino_Run2015C_PromptReco_50ns_DoubleEG.root");
    myProject->AddDataFile(pathData + "50ns/latino_Run2015C_PromptReco_50ns_MuonEG.root");
    myProject->AddDataFile(pathData + "50ns/latino_Run2015C_PromptReco_50ns_SingleElectron.root");
    myProject->AddDataFile(pathData + "50ns/latino_Run2015C_PromptReco_50ns_SingleMuon.root");
    myProject->AddDataFile(pathData + "50ns/latino_Run2015C_PromptReco_50ns_DoubleMuon.root");
    */
    myProject->AddDataFile(pathData + "50ns/latino_DoubleEG.root");
    myProject->AddDataFile(pathData + "50ns/latino_MuonEG.root");
    myProject->AddDataFile(pathData + "50ns/latino_SingleElectron.root");
    myProject->AddDataFile(pathData + "50ns/latino_SingleMuon.root");
    myProject->AddDataFile(pathData + "50ns/latino_DoubleMuon.root");

    isdata             = true;
    nEventsInTheSample = 128512; 
    xSection           = 12.461;
    whichRun           = 2; 
  }

  else if (data=="WW50") {

    myProject->AddDataFile(path50 + "latino_WWTo2L2Nu_NLL.root");

    isdata             = false;
    nEventsInTheSample = 376031;
    xSection           = 10.481;
    whichRun           = 2; 
  }

  else if (data=="WJets50") {

    myProject->AddDataFile(path50 + "latino_WJetsToLNu.root");

    isdata             = false;
    nEventsInTheSample = 10249764;
    xSection           = 20508.9;
    whichRun           = 2;
  }

  else if (data=="VV50") {
  
    myProject->AddDataFile(path50 + "latino_WZ.root");
    myProject->AddDataFile(path50 + "latino_ZZ.root");

    isdata             = false;
    nEventsInTheSample = 266293 + 204666;
    xSection           = 31.8 + 66.1;
    whichRun           = 2;
  }
  
  else if (data == "TTJets50") {

    myProject->AddDataFile(path50 + "latino_TTJets.root");

    isdata             = false;
    nEventsInTheSample = 2562273;
    xSection           = 831.76;
    whichRun           = 2;
  }
  /* GEN_weight_SM still not available
  else if (data == "SingleTop50") {

    myProject->AddDataFile(path50 + "latino_ST_t-channel.root");

    isdata             = false;
    nEventsInTheSample = 2309030;
    xSection           = 87.315;
    whichRun           = 2;
  }
  */
  else if (data=="Top50") {

    //myProject->AddDataFile(path50 + "latino_TTJets.root");
    //myProject->AddDataFile(path50 + "latino_ST_t-channel.root");
    myProject->AddDataFile(path25 + "latino_ST_t-channel_top.root");     //25ns sample!!!
    myProject->AddDataFile(path25 + "latino_ST_t-channel_antitop.root"); //25ns sample!!!
    myProject->AddDataFile(path50 + "latino_ST_tW_antitop.root");
    myProject->AddDataFile(path50 + "latino_ST_tW_top.root");

    isdata             = false;
    nEventsInTheSample = 2562273 * 3;
    xSection           = 70.69 + 35.6 + 35.6;
    whichRun           = 2;
  }

  else if (data=="DY50") {
    
    myProject->AddDataFile(path50 + "latino_DYJetsToLL_M-50.root");

    isdata             = false;
    nEventsInTheSample = 11934836;
    xSection           = 6025.2;
    whichRun           = 2;
  }
  
  else if (data=="HWW50") {

    myProject->AddDataFile(path25 + "latino_GluGluHToWWTo2L2Nu_M125.root");

    isdata             = false;
    nEventsInTheSample = 475039 + 473153;
    xSection           = 71.6;
    whichRun           = 2;
  }

  else{
    cout<<"************************************************************"<<endl;
    cout<<"I can't find the sample you are asking for. Please try again"<<endl;
    cout<<"************************************************************"<<endl;

    return;
  }

  //Number of events to process
  myProject->SetNEvents(nEventsToProcess);

  ///////////////////////////////
  // INPUT PARAMETERS
 
  myProject->SetInputParam("IsDATA",        isdata);
  myProject->SetInputParam("Signal",        data);
  myProject->SetInputParam("XSection",      xSection);
  myProject->SetInputParam("Luminosity",    luminosity);
  myProject->SetInputParam("NEvents",       nEventsInTheSample);
  myProject->SetInputParam("luminosityPU",  19468.3);  
  myProject->SetInputParam("WhichRun",      whichRun);
  myProject->SetInputParam("Debug",         debug);
  myProject->SetInputParam("Report",        doReport);
  myProject->SetInputParam("FlavorChannel", flavorChannel);
  myProject->SetInputParam("SameSign",      sameSign);
  myProject->SetInputParam("OutPath",       outPath);
  myProject->SetInputParam("MuonID",        muonID);

  ///////////////////////////////
  // OUTPUT FILE NAME
  // Specify the name of the file where you want your histograms to be saved

  if (data.Contains("25")){
    TString name = data;
    name.Remove(name.Length() - 2, 2);
    myProject->SetOutputFile(outPath + "25ns/" + name + ".root");
  }

  else if (data.Contains("50")){
    TString name = data;
    name.Remove(name.Length() - 2, 2);
    myProject->SetOutputFile(outPath + "50ns/" + name + ".root");
  }

  ///////////////////////////////
  // SELECTOR AND PACKAGES

  //Add the core selector, this one is essential 
  myProject->AddSelectorPackage("WWAnalysisSelector");

  //Add the additional selectors if needed (and if available!)

  // RUN THE ANALYSIS
  // ================

  myProject->Run();

  /////////////////////////////////////////////////////////////////////////

}
