//root -l -q 'yields.C("MuMu")'

#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

std::ofstream inFile("yields.txt",std::ios::out);

//defining processes (signal & backgrounds)
const int nProcesses = 3;

enum {iWW, iTT, iWJets};

TFile *input[nProcesses];
TH1F  *histo[nProcesses];

TString process[nProcesses];

process[iWW]     = "WW50";
process[iTT]     = "TTbar50";
process[iWJets]  = "WJets50";

Float_t yield[nProcesses];

const int nCuts = 10;

//defining cut levels
TString cutLevel[nCuts];
cutLevel[0] = "Two Leptons";
cutLevel[1] = "No Extra Leptons";
cutLevel[2] = "pfMet > 20~GeV";
cutLevel[3] = "m$_{ll}$ > 12~GeV";
cutLevel[4] = "|m$_{ll}$ - m$_Z$| > 15~GeV";
cutLevel[5] = "mpMet > 20~GeV";
cutLevel[6] = "$\\Delta_{\phi}$ veto";
cutLevel[7] = "p$_T^{ll}$ > 30~GeV (in OF p$_T^{ll}$ > 45~GeV)";
cutLevel[8] = "B-veto";
cutLevel[9] = "Ht < 250~GeV";

TString cutHisto[nCuts];
cutHisto[0] = "hDeltaPhiLeptonsTwoLeptonsLevel";
cutHisto[1] = "hWExtraLepton";
cutHisto[2] = "hWMetCut";
cutHisto[3] = "hWLowMinv";
cutHisto[4] = "hWZVeto";
cutHisto[5] = "hWpMetCut";
cutHisto[6] = "hWDeltaPhiJet";
cutHisto[7] = "hWPtll";
cutHisto[8] = "hWTopTagging";
cutHisto[9] = "hDeltaPhiLeptonsWWLevel3";

void yields(TString flavourChannel = ""){

  if (flavourChannel != "All" && flavourChannel !=  "OF" && flavourChannel !=  "SF" && flavourChannel !=  "MuMu" && flavourChannel !=  "EE" && flavourChannel !=  "EMu" && flavourChannel !=  "MuE"){

    cout<<"***********************************************************************************************************"<<endl;
    cout<<"Please select a valid flavour channel to look at: 'All' or 'OF' or 'SF' or 'MuMu' or 'EE' or 'EMu' or 'MuE'"<<endl;
    cout<<"For example: root -l -q 'yields.C(\"OF\")'"<<endl;
    cout<<"***********************************************************************************************************"<<endl;
    return;
  }

  TString path = "rootFiles/" + flavourChannel + "/";
  
  inFile<<"\\begin{tabular}{cSSSS}"<<endl;
  cout<<"\\begin{tabular}{cSSSS}"<<endl;
  inFile<<"\\toprule"<<endl;
  cout<<"\\toprule"<<endl;
  inFile<<"Cut Level";
  cout<<"Cut Level";
  for (int q = 0; q < nProcesses; ++q){
    inFile<<" & {"<<process[q]<<"}";
    cout<<" & {"<<process[q]<<"}";
  }			
  inFile<<" \\\\ "<<endl;	
  cout<<" \\\\ "<<endl;
  inFile<<"\\midrule"<<endl;
  cout<<"\\midrule"<<endl;

  for(int c = 0; c < nCuts; ++c){
    
    //signal and background files and histograms
    for (int ip = 0; ip < nProcesses; ++ip){
      TString fileName = path + process[ip] + ".root";
      input[ip] = new TFile(fileName, "read");
      
      if (input[ip] -> GetListOfKeys() -> Contains(cutHisto[c]) == 0){
	cout<<"I cannot find the histogram named "<<cutHisto[c]<<". I'll skip it."<<endl;
	return 0;
      }
      
      histo[ip]  = (TH1F*) input[ip] -> Get(cutHisto[c]);
    }
    
    //begin the table
    inFile<<cutLevel[c];
    cout<<cutLevel[c];
    //number of events
    for(int i = 0; i < nProcesses; ++i){
      yield[i] = TMath::Nint(histo[i] -> Integral()); 
      inFile<<" & "<<yield[i];
      cout<<" & "<<yield[i];
    }
    inFile<<"\\\\"<<endl;
    cout<<"\\\\"<<endl;
  }
  
  //end of the table
  inFile<<"\\bottomrule"<<endl;
  cout<<"\\bottomrule"<<endl;
  //inFile<<"Percentage"<<" & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;
  inFile<<"\\end{tabular}"<<endl;
  cout<<"\\end{tabular}"<<endl;
  inFile.close();

  gSystem -> Exec("mv yields.txt distributions/" + flavourChannel );
}
