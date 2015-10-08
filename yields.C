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
#include <iomanip>
#include <vector>

std::ofstream inFile("yields.txt",std::ios::out);

std::vector<float> sigVector;
std::vector<float> dataVector;
std::vector<float> totMC;

const int nProcesses = 8;
TString process[nProcesses];

enum {iWW, iVV, iTop, iTTJets, iDY, iWJets, iHWW, iData};

process[iWW]     = "WW";
process[iWJets]  = "WJets";
process[iVV]     = "VV";
process[iTop]    = "Top";
process[iTTJets] = "TTJets";
process[iDY]     = "DY";
process[iHWW]    = "HWW";
process[iData]   = "Data2015";

/*
const int nProcesses = 6;
TString process[nProcesses];

enum {iWW, iVV, iTop, iDY, iWJets, iHWW};

process[iWW]    = "WW";
process[iVV]    = "VV";
process[iTop]   = "Top";
process[iDY]    = "DY";
process[iWJets] = "WJets";
process[iHWW]   = "HWW";
*/
TFile *input[nProcesses];
TH1F  *histo[nProcesses];

Float_t yield[nProcesses];

const int nCuts = 15;

//TGraphErrors *gSignificance = new TGraphErrors();
TH1F *gSignificance = new TH1F("Significance","Significance",20,0,20);

//defining cut levels
TString cutLevel[nCuts];
cutLevel[0]  = "Two Leptons";
cutLevel[1]  = "No Extra Leptons";
cutLevel[2]  = "pfMet > 20~GeV";
cutLevel[3]  = "m$_{ll}$ > 12~GeV";
cutLevel[4]  = "|m$_{ll}$ - m$_Z$| > 15~GeV and MetVar > 45";
cutLevel[5]  = "mpMet > 20~GeV";
cutLevel[6]  = "$\\Delta_{\phi}$ veto";
cutLevel[7]  = "p$_T^{ll}$ >";
cutLevel[8]  = "B-veto [IP and tche]";
cutLevel[9]  = "B-veto [New IP]";
cutLevel[10] = "B-veto [csvv2ivf Loose]";
cutLevel[11] = "B-veto [csvv2ivf Medium]";
cutLevel[12] = "B-veto [csvv2ivf Tight]";
cutLevel[13] = "B-veto [csvv2ivf Loose and Mu]";
cutLevel[14] = "Ht < 237~GeV";
             
TString cutHisto[nCuts];
cutHisto[0]  = "hWTrigger";
cutHisto[1]  = "hWExtraLepton";
cutHisto[2]  = "hWMetCut";
cutHisto[3]  = "hWLowMinv";
cutHisto[4]  = "hWZVeto";
cutHisto[5]  = "hWpMetCut";
cutHisto[6]  = "hWDeltaPhiJet";
cutHisto[7]  = "hWPtll";
cutHisto[8]  = "hbVetoOld3";
cutHisto[9]  = "hbVetoMu3";
cutHisto[10] = "hbVetoCsvv2ivfLoose3";
cutHisto[11] = "hbVetoCsvv2ivfMedium3";
cutHisto[12] = "hbVetoCsvv2ivfTight3";
cutHisto[13] = "hbVetoCsvv2ivfLooseAndMu3";
cutHisto[14] = "hDeltaPhiLeptonsWWLevel3";
             
void yields(TString flavourChannel = "",
	    TString muonID         = "",
	    TString bunch          = ""
	    ){
  
  if( muonID != "MediumID" && muonID != "TightID" && muonID != "TightIDTighterIP")  
    if (flavourChannel != "All" && flavourChannel !=  "OF" && flavourChannel !=  "SF" && flavourChannel !=  "MuMu" && flavourChannel !=  "EE" && flavourChannel !=  "EMu" && flavourChannel !=  "MuE" && bunch != "50ns" && bunch != "25ns"){
      
      cout<<"**************************************************************************************************************"<<endl;
      cout<<"Please select a valid flavour channel to look at: 'All' or 'OF' or 'SF' or 'MuMu' or 'EE' or 'EMu' or 'MuE'..."<<endl;
      cout<<"... a valid muonID: 'MediumID', 'TightID' or 'TightIDTighterIP'..."<<endl;
      cout<<"... and which bunch spacing you want to analyze: '25ns' or '50ns'"<<endl;
      cout<<"For example: root -l -q 'yields.C(\"OF\",\"MediumIDTighterIP\",\"50ns\")'"<<endl;
      cout<<"**************************************************************************************************************"<<endl;
      return;
    }
  
  if (flavourChannel == "OF" || flavourChannel == "EMu" || flavourChannel == "MuE")
    cutLevel[7] = cutLevel[7] + " 30~GeV";
  
  else if (flavourChannel == "SF" || flavourChannel == "MuMu" || flavourChannel == "EE")
    cutLevel[7] = cutLevel[7] + " 45~GeV";
  
  else if (flavourChannel == "All"){
    cutLevel[6] = cutLevel[6] + " (SF only)";
    cutLevel[7] = "p$_T^{ll}$ > 30~GeV (45~GeV) in OF (SF)";
  }
  
  TString path = "rootFiles/" + flavourChannel + "/" + muonID + "/" + bunch + "/";
  
  inFile<<"\\begin{tabular}{cSSSSSSSS}"<<endl;
  cout<<"\\begin{tabular}{cSSSSSSSS}"<<endl;
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
    if (flavourChannel == "OF" || flavourChannel == "EMu" || flavourChannel == "MuE") 
      if (c == 6 || c == 4) {
	sigVector.push_back(0);
	dataVector.push_back(0);
	totMC.push_back(0);
	continue;
      }

    inFile<<cutLevel[c];
    cout<<cutLevel[c];
    
    Float_t den = 0.;
    
    //number of events
    for(int i = 0; i < nProcesses; ++i){
      
      yield[i] = /*TMath::Nint*/(histo[i] -> Integral()); 
      if (i != iData)
	den += yield[i];
      else dataVector.push_back(yield[i]);
      std::cout << std::fixed;
    
      inFile<<" & "<<setprecision(2)<<yield[i];
      cout<<"   & "<<setprecision(2)<<yield[i];
    }
    gSignificance -> SetBinContent(c,yield[iWW] / sqrt(den));
    sigVector.push_back(yield[iWW] / sqrt(den));
    totMC.push_back(den);
    TString dump = cutLevel[c];
    dump.Remove(0,8);
    dump.Remove(dump.Length()-1,1);
    gSignificance -> GetXaxis() -> SetBinLabel(c, dump);
    inFile<<"\\\\"<<endl;
    cout<<"\\\\"<<endl;
  }
  
  //end of the table
  inFile<<"\\bottomrule"<<endl;
  cout<<"\\bottomrule"<<endl;
  //inFile<<"Percentage"<<" & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;
  inFile<<"\\end{tabular}"<<endl;
  cout<<"\\end{tabular}"<<endl;

  //Significance Tabular
  inFile<<"\\begin{center}"<<endl;
  cout<<"\\begin{center}"<<endl;
  inFile<<"\\begin{tabular}{cSSSSSSS}"<<endl;
  cout<<"\\begin{tabular}{cSSSSSSS}"<<endl;
  inFile<<"\\toprule"<<endl;
  cout<<"\\toprule"<<endl;
  inFile<<"Cut Level";
  cout<<"Cut Level";
  inFile<<" & {Significance}";
  cout<<" & {Significance}";
  inFile<<" & {Data2015}";
  cout<<" & {Data2015}";
  inFile<<" & {Total MC}";
  cout<<" & {Total MC}";
  inFile<<" \\\\ "<<endl;	
  cout<<" \\\\ "<<endl;
  inFile<<"\\midrule"<<endl;
  cout<<"\\midrule"<<endl;
  for(int i = 0; i < nCuts; ++i){
    if (flavourChannel == "OF" || flavourChannel == "EMu" || flavourChannel == "MuE")
      if (i == 6 || i == 4) continue;
    
    cout<<cutLevel[i];
    inFile<<cutLevel[i];
    cout<<" & "<<sigVector.at(i);
    inFile<<" & "<<sigVector.at(i);
    cout<<" & "<<dataVector.at(i);
    inFile<<" & "<<dataVector.at(i);
    cout<<" & "<<totMC.at(i);
    inFile<<" & "<<totMC.at(i);
    inFile<<"\\\\"<<endl;
    cout<<"\\\\"<<endl;
  }
  inFile<<"\\bottomrule"<<endl;
  cout<<"\\bottomrule"<<endl;
  inFile<<"\\end{tabular}"<<endl;
  cout<<"\\end{tabular}"<<endl;
  inFile<<"\\end{center}"<<endl;
  cout<<"\\end{center}"<<endl;
  inFile.close();


  //Significance Graph Cosmetics
  gSignificance -> GetXaxis() -> SetNdivisions(-414);
  gSignificance -> GetXaxis() -> SetLabelSize(0.04);
  gSignificance -> GetYaxis() -> SetLabelSize(0.04);
  gSignificance -> GetXaxis() -> SetRangeUser(7,13);
  gSignificance -> GetYaxis() -> SetRangeUser(0.,3.5);
  gSignificance -> GetXaxis() -> CenterLabels();
  gSignificance -> GetXaxis() -> LabelsOption("h");
  gSignificance -> SetMarkerStyle(8);
  gSignificance -> SetStats(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1 -> cd();
  gSignificance->Draw("p");
  c1 -> Print("Bveto.pdf","pdf");

  gSystem -> Exec("mv yields.txt distributions/" + bunch + "/" + flavourChannel + "/" + muonID);
  gSystem -> Exec("mv Bveto.pdf  distributions/" + bunch + "/" + flavourChannel + "/" + muonID);
}
