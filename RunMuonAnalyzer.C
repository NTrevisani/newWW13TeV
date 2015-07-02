///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		     Float_t luminosity  //In pb-1
		     ) {
 
  gSystem->Load("libPAF.so");   

  TString path25      = "/gpfs/csic_projects/cms/trevisanin/newLatino/25ns/";
  TString path50      = "/gpfs/csic_projects/cms/trevisanin/newLatino/50ns/";
  TString pathOld     = "/gpfs/csic_projects/cms/trevisanin/newLatino/";

  TString outPath     = "rootFiles/" + flavorChannel + "/";

  gSystem -> Exec("mkdir rootFiles");
  gSystem -> Exec("mkdir rootFiles/" + flavorChannel);

  // Manual Input Parameters
  bool     debug            = false;   //For verbose while debugging
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

  else {
    cout<<"Please select a valid PAF operating mode: Sequential or Cluster"<<endl;
    return; 
  }

  // Set default name and subdirectory of Trees
  myProject->SetDefaultTreeName("latino");

  ////////////////////////////////////////////////
  // ADD SAMPLES

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
    nEventsInTheSample = 132180; 
    xSection           = 20508.9;
    whichRun           = 2; 
  }

  else if (data=="WW50") {

    myProject->AddDataFile(path50 + "latino_WWTo2L2Nu.root");

    isdata             = false;
    nEventsInTheSample = 128512;
    xSection           = 12.461;
    whichRun           = 2;
  }

  else if (data=="WJets50") {

    myProject->AddDataFile(path50 + "latino_WJetsToLNu.root");

    isdata             = false;
    nEventsInTheSample = 132180;
    xSection           = 20508.9;
    whichRun           = 2;
  }

  else if (data=="TTbar50") {

    myProject->AddDataFile(path50 + "latino_TTTo2L2Nu.root");

    isdata             = false;
    nEventsInTheSample = 2309030;
    xSection           = 87.315;
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

  ///////////////////////////////
  // OUTPUT FILE NAME
  // Specify the name of the file where you want your histograms to be saved

  myProject->SetOutputFile(outPath + data + ".root");

  ///////////////////////////////
  // SELECTOR AND PACKAGES

  //Add the core selector, this one is essential 
  myProject->AddSelectorPackage("CoreMuonSelector");

  //Add the additional selectors if needed (and if available!)

  // RUN THE ANALYSIS
  // ================

  myProject->Run();

  /////////////////////////////////////////////////////////////////////////

}
