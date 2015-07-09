/*
root -l -b -q 'ROC.C("png","WW50")'
*/

#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "math.h"
#include "fstream.h"
#include <vector>
#include <string>

const int nBkg = 2;
const int nSig = 1;

//signal files
TFile *fSig[nSig]; 

//defining signal names
TString sigName[nSig];
sigName[0] = "WW50";

//defining signal histogram
TH1F *hSig[nSig];

//defining background names
TString bkgName[nBkg];
bkgName[0] = "TTbar50"; 
bkgName[1] = "WJets50";  

//background files
TFile *fBkg[nBkg];

//defining background histograms
TH1F *hBkg[nBkg];

//defining graphs
TGraph *g  = new TGraph();
TGraph *g2 = new TGraph();
TGraph *g3 = new TGraph();
TGraph *g4 = new TGraph();

//file where I'll put the cut values for each variable
std::ofstream cutFile("info.txt",std::ios::out);

std::vector<float> significances;
std::vector<float> cuts;

typedef std::vector<std::string> str_vec_t;
//std::vector<std::string> variables;
//std::vector<std::string> direction;

str_vec_t variables;
str_vec_t direction;

float GetXofTheMax(TGraph* grafico, int endPoint){
  float value = 0.;
  float xmax = 0.;
  for (int  i = 0; i < endPoint; ++i){
    if (grafico->Eval(i) > value){
      value = grafico->Eval(i);
      xmax = i;
    }
  }
  return xmax;
}

float GetMax(TGraph* grafico, int endPoint){
  float value = 0.;
  for (int  i = 0; i < endPoint; ++i){
    if (grafico->Eval(i) > value){
      value = grafico->Eval(i);
    }
  }
  return value;
}

float FitBestCut(TGraph* grafico, int end){
  float mean = GetXofTheMax(grafico,end);
  TF1* curva = new TF1("curva","pol5",0.,3000.);
  grafico->Fit("curva","","",mean/2,1.5*mean);
  float max = -9999;
  float maxX = -9999;
  float passo = mean/10000.;
  for (int i = 0; i < 10000; ++i){
    if(curva->Eval(mean/2 + passo*i) > max){
      max = curva->Eval(mean/2 + passo*i);
      maxX = mean/2 + passo*i;
    }
  }
  return maxX;
}

float FitBestY(TGraph* grafico, int end){
  float mean = GetXofTheMax(grafico,end);
  TF1* curva = new TF1("curva","pol5",0.,3000.);
  grafico->Fit("curva","","",mean/2,1.5*mean);
  float max = -9999;
  float passo = mean/10000.;
  for (int i = 0; i < 10000; ++i){
    if(curva->Eval(mean/2 + passo*i) > max){
      max = curva->Eval(mean/2 + passo*i);
    }
  }
  return max;
}

int doTheROC(TString variable = "", 
	     TString typeOfCut = "",
	     TString max_ = "",
	     TString printFormat = "",
	     TString text = ""
	     ){

  gSystem->Exec("mkdir ROCfolder");
  gSystem->Exec("mkdir ROCfolder/" + printFormat);

  if(typeOfCut != ">" && typeOfCut != "<"){
    cout<<"Please tell me if you want to keep the events above the optimal value I'll find ('>') or below it ('<')"<<endl;
    return 0;
  }

  if(printFormat != "png" && printFormat != "pdf" && printFormat != "C"){
    cout<<"Plese select a valid output format for the plots: pdf, png or C"<<endl;
    return 0;
  }

  Float_t max = 0;
  if (max_.IsDigit())
    max = max_.Atof();

  //signal files and histograms
  for(int j = 0; j < nSig; ++j){
    fSig[j] = new TFile("rootFiles/OF/" + sigName[j]  + ".root","read");
    hSig[j] = (TH1F*)fSig[j] -> Get(variable);
  }

  //background files and histograms
  for(int i = 0; i < nBkg; ++i){
    fBkg[i] = new TFile("rootFiles/OF/" + bkgName[i]  +".root","read");
    hBkg[i]  = (TH1F*)fBkg[i] -> Get(variable);
  }

  //total number of signal events (no cuts applied)
  Float_t totSig = hSig[0]->Integral();

  //total number of background events (no cuts applied)
  Float_t totBkg = 0;
  for(int q = 0; q < nBkg; ++q)
    totBkg = totBkg + hBkg[q] ->Integral();

  //defining variables
  Float_t effSig = 0.;
  Float_t effBkg = 0.;

  //main loop: scanning efficiencies
  for (Int_t i = 0; i < 3000; ++i){
    
    //calculating signal efficiency
    Float_t sum = 0.;
    if (typeOfCut == "<")
      sum = hSig[0]->Integral(0.,i);                              
    else if (typeOfCut == ">")
      sum = hSig[0]->Integral(3000. - i, 3000.);                  
    if(totSig != 0)
      effSig = sum / totSig;
      
    //calculating background efficiency (and rejection)
    Float_t bkgInt = 0;
    for(int pp = 0; pp < nBkg; ++pp)
      if (typeOfCut == "<")
	bkgInt = bkgInt + hBkg[pp]->Integral(0.,i);            
      else if (typeOfCut == ">")      
	bkgInt = bkgInt + hBkg[pp]->Integral(3000. - i, 3000.);

    if(totBkg != 0)
      effBkg = bkgInt / totBkg; 

    //filling graphs
    g ->SetPoint(i, effSig, 1 - effBkg);

    if (typeOfCut == "<"){
      Float_t den = (bkgInt + hSig[0]->Integral(0.,i));           
      if (den != 0)
	g2->SetPoint(i, i, (hSig[0]->Integral(0.,i)) / sqrt(den));
      g3->SetPoint(i, i, effSig);                              
      g4->SetPoint(i, i, 1 - effBkg); 
    }   

    else if (typeOfCut == ">"){
      Float_t den = (bkgInt + hSig[0]->Integral(3000. - i, 3000.));   
      if (den != 0)
	g2->SetPoint(i, 3000. - i, (hSig[0]->Integral(3000. - i, 3000.)) / sqrt(den));
      g3->SetPoint(i, 3000. - i, effSig);
      g4->SetPoint(i, 3000. - i, 1. - effBkg);
    }
  }

  TCanvas *cb = new TCanvas("cb","cb",600.,600.);
  cb -> cd();

  TPad *padb = new TPad("padb","padb",0.,0.,1.,1.);
  padb->SetLeftMargin(0.15);
  padb->SetRightMargin(0.15);
  padb->SetBottomMargin(0.15);
  padb->Draw();
  padb->cd();

  hSig[0] -> Rebin(10);
  hSig[0] -> Scale (1. / hSig[0] -> Integral());
  hSig[0] -> GetYaxis()->SetRangeUser(0.,0.0180);
  hSig[0] -> SetLineWidth(5);
  hSig[0] -> SetStats(0);

  hSig[0] -> GetXaxis()->SetTitle(text);
  hSig[0] -> SetTitle("Normalized Distributions for the Sum of Leptons p_{T}");
  hSig[0]->GetXaxis()->SetTitleSize(0.05);
  hSig[0]->GetXaxis()->SetTitleOffset(1.2);
  hSig[0]->GetXaxis()->SetLabelSize(0.05);
  hSig[0]->GetYaxis()->SetTitleOffset(1.4);
  hSig[0]->GetYaxis()->SetTitleSize(0.05);
  hSig[0]->GetYaxis()->SetLabelSize(0.05);
  hSig[0]->GetXaxis()->SetNdivisions(408);

  hSig[0] -> Draw();

  TLegend* legb = new TLegend(0.60,0.55,0.75,0.85);
  legb->SetTextSize(0.04);
  legb->SetFillColor(kWhite);
  legb->SetLineColor(kWhite);

  legb -> AddEntry(hSig[0],"signal","l");

  for(int b = 0; b < nBkg; ++b){
    hBkg[b] -> Rebin(10);
    hBkg[b] -> SetStats(0);
    hBkg[b] -> SetLineWidth(5);
    hBkg[b] -> Scale (1. / hBkg[b] -> Integral());
    hBkg[b] -> GetYaxis()->SetRangeUser(0.,0.018);
    hBkg[b] -> SetLineColor(b + 1);
    legb -> AddEntry(hBkg[b], bkgName[b],"l");
    cout<<hBkg[b]->GetTitle()<<endl;
    hBkg[b] -> Draw("same");
  }
  legb->Draw();

  g->SetTitle("ROC");
  g->GetXaxis()->SetTitle("Signal Efficiency");
  g->GetYaxis()->SetTitle("Background Rejection");
  g->SetLineWidth(5);  
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(1.2);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetXaxis()->SetRangeUser(0.,1.);
  g->GetYaxis()->SetRangeUser(0.,1.);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);

  TGraph *p = new TGraph();
  p->SetPoint(0.,0.84,0.55);
  p->SetMarkerColor(kRed);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(2);
  
  g2->SetTitle("Significance");
  g2->GetXaxis()->SetTitle(text);
  g2->SetLineWidth(5);
  g2->GetYaxis()->SetRangeUser(0.,1.5*GetMax(g2,max));//0.018);
  g2->GetXaxis()->SetRangeUser(0.,max);
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetXaxis()->SetTitleOffset(1.2);
  g2->GetXaxis()->SetLabelSize(0.05);
  g2->GetYaxis()->SetTitle("S / #sqrt{S + B}");
  g2->GetYaxis()->SetTitleOffset(2.0);
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->SetLabelSize(0.05);

  //  Float_t drawCut = FitBestCut(g2,max);
  Float_t drawCut = GetXofTheMax(g2,max);
  cout<<"I'd like to cut here: "<<TMath::Nint(drawCut)<<endl;
  //Float_t significanceCut = FitBestY(g2,max);
  Float_t significanceCut = g2->Eval(TMath::Nint(drawCut));
  cout<<"Doing so the Significance would be "<<significanceCut<<endl;
  //cutFile<<variable<<" "<<typeOfCut<<" "<<TMath::Nint(drawCut)<<" -> Significance  "<<significanceCut<<endl;

  std::string typeOfCutStd = typeOfCut;
  std::string variableStd = variable;

  cuts.push_back(drawCut);
  significances.push_back(significanceCut);
  variables.push_back(variableStd);
  direction.push_back(typeOfCutStd);  

  TLine *l = new TLine(drawCut,0.,drawCut,1.5*GetMax(g2,max));
  l->SetLineColor(kViolet);
  l->SetLineWidth(5);

  g3->SetLineColor(kRed);
  g3->SetLineWidth(5);
  g3->SetTitle("Signal Efficiency and Background Rejection");
  g3->GetXaxis()->SetTitle(text);
  g3->GetYaxis()->SetTitle("Signal Efficiency");
  g3->GetXaxis()->SetRangeUser(0.,max);  
  g3->GetYaxis()->SetRangeUser(0.,1.);
  g3->GetXaxis()->SetLabelOffset(0.02);
  g3->GetXaxis()->SetTitleSize(0.05);
  g3->GetXaxis()->SetTitleOffset(1.3);
  g3->GetXaxis()->SetLabelSize(0.05);
  g3->GetYaxis()->SetTitleOffset(1.4);
  g3->GetYaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->SetLabelSize(0.05);
  g4->SetLineColor(kBlack);
  g4->SetLineWidth(5);

  TLegend* leg = new TLegend(0.60,0.50,0.75,0.70);
  leg->AddEntry(g3,"Sig Eff","l");
  leg->AddEntry(g4,"Bkg Rej","l");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);

  TLine *ll = new TLine(drawCut,0.,drawCut,1.);
  ll->SetLineColor(kBlue);
  ll->SetLineWidth(5);

  TCanvas *c1 = new TCanvas("ROC","ROC",600.,600.);
  c1->cd();

  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.15);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  g->Draw("AL");
  //p->Draw("Psame");

  c1->Print("RocAll" + variable + "." + printFormat,printFormat); 
  gSystem->Exec("mv RocAll" + variable + "." + printFormat + " ROCfolder/" + printFormat);
  
  TCanvas *c2 = new TCanvas("Significance","Significance",600.,600.);
  c2->cd();

  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1.0, 1.0);
  pad2->SetLeftMargin(0.20);
  pad2->SetBottomMargin(0.15);
  pad2->Draw();
  pad2->cd();

  g2->Draw("AL");
  l->Draw("same");
  //  g4->Draw("same");

  c2->Print("SignificanceAll" + variable + "." + printFormat, printFormat);  
  gSystem->Exec("mv SignificanceAll" + variable + "." + printFormat + " ROCfolder/" + printFormat);

  TCanvas *c3 = new TCanvas("Efficiencies","Efficiencies",600.,600.);
  c3->cd();

  TPad* pad3 = new TPad("pad3", "pad3", 0., 0., 1.0, 1.0);
  pad3->SetLeftMargin(0.15);
  pad3->SetRightMargin(0.15);
  pad3->SetBottomMargin(0.15);
  pad3->Draw();
  pad3->cd();

  g3->Draw("AL");
  g4->Draw("Lsame");
  ll->Draw("same");
  pad3->Update();
  TGaxis *axis = new TGaxis(pad3->GetUxmax(),pad3->GetUymin(),pad3->GetUxmax(),pad3->GetUymax(),0,1,510,"+L");
  axis->SetLabelSize(g3->GetXaxis()->GetLabelSize());
  axis->SetLabelFont(g3->GetXaxis()->GetLabelFont());
  axis->SetTitleSize(g3->GetXaxis()->GetTitleSize());
  axis->SetTitleFont(g3->GetXaxis()->GetTitleFont());
  axis->SetTitleOffset(1.4);
  axis->SetTitle("Background Rejection");
  axis->Draw("same");
  leg->Draw("same");

  c3->Print("EfficienciesAll" + variable + "." + printFormat,printFormat);
  gSystem->Exec("mv EfficienciesAll" + variable + "." + printFormat + " ROCfolder/" + printFormat);
  
  TCanvas *c4 = new TCanvas("AllTogether","AllTogether",1800.,600.);
  c4->cd();

  c4->Divide(3,1);

  c4->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  g->Draw();
  c4->cd(2);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  g2->Draw();
  l->Draw("same");
  c4->cd(3);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  g3->Draw();
  g4->Draw("Lsame");
  axis->Draw("same");
  leg->Draw("same");
  ll->Draw("same");

  c4->Print(variable + "." + printFormat,printFormat);
  gSystem->Exec("mv " + variable + "." + printFormat + " ROCfolder/" + printFormat);

  delete c1;
  delete c2;
  delete c3;
  g ->Clear();
  g2->Clear();
  g3->Clear();
  g4->Clear();
  /*  for(int d = 0; d < nSig; ++d){
    delete fSig[d];
    delete hSig[d];
  }
  for(int e = 0; e < nBkg; ++e){
    delete fBkg[d];
    delete hBkg[d];
  }
  */
  return 0;
}

void ROC(TString printMode = ""){

  //  if (variable == "" || typeOfCut == "" || max_ == "" || printFormat == "" || text == ""){
  if(printMode == ""){
    cout<<"********************************************************************************************"<<endl;
    cout<<"Please select a valid set of input options:"<<endl;
    //cout<<"variable: name of the histogram you want to analyze"<<endl;
    //cout<<"typeOfCut: select the part of the histogram you want to KEEP, '>' or '<' wrt the cut applied"<<endl;
    //cout<<"max_: the maximum value you want the code to scan the histogram (no more than 3000)"<<endl;
    cout<<"printFormat: 'C', 'png' or 'pdf'"<<endl;
    //cout<<"text: the x-axis tag (usually the name of the variable and the units)"<<endl;
    cout<<"Remember that some pieces of information have to be added to the 'letsROC.txt file:'"<<endl;
    cout<<"'variable', 'typeOfCut', 'max_' and 'text'"<<endl;
    cout<<"Launch the code as: root -l -b -q 'ROC.C('printFormat')'"<<endl;
    cout<<"For example: root -l -b -q 'ROC.C(\"pdf\")'"<<endl;
    cout<<"********************************************************************************************"<<endl;
    return 0;
  }

  TString var = ""; 
  TString cutType = "";
  TString rightBound = "";
  TString tagText = "";

  //input file
  ifstream inFile("letsROC.txt");
  std::string line;
  
  while (getline(inFile,line)){
    std::ofstream outFile("out.tmp",std::ios::out);
    outFile<<line<<endl;
    outFile.close();
    std::ifstream input;
    input.open("out.tmp",std::ios::in);
    input >> var >> cutType >> rightBound >> tagText;
    input.close();
    //    cout<< var <<" "<<" "<< cutType <<" "<< rightBound <<" "<<printMode<<" "<< tagText <<endl;
    doTheROC(var, cutType, rightBound, printMode, tagText);
  } 

  float helpSig;
  std::string helpVar;
  std::string helpDir;
  float helpCut;

  //order the variables in significance
  for (int o = 0; o < significances.size(); ++o)
    for (int p = o; p < significances.size(); ++p)
      if(significances.at(o) < significances.at(p)){
	helpSig = significances.at(o);
	significances.at(o) = significances.at(p);
	significances.at(p) = helpSig;
	helpVar = variables.at(o);
	variables.at(o) = variables.at(p);
	variables.at(p) = helpVar;
	helpDir = direction.at(o);
	direction.at(o) = direction.at(p);
	direction.at(p) = helpDir;
	helpCut = cuts.at(o);
	cuts.at(o) = cuts.at(p);
	cuts.at(p) = helpCut;
      }
  for (int j = 0; j < significances.size(); ++j)
    cutFile<<variables.at(j)<<" "<<direction.at(j)<<" "<<TMath::Nint(cuts.at(j))<<" -> Significance  "<<significances.at(j)<<endl;

  inFile.close();
  cutFile.close();
  gSystem->Exec("rm out.tmp");
  gSystem->Exec("mv info.txt ROCfolder/" + printMode);
}
