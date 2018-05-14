#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <THStack.h>
#include <TList.h>
#include <TFrame.h>
#include <TPostScript.h>
#include <TPaveLabel.h>
#include <TText.h>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TColor.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;
#endif

int main(int argc, const char* argv[]){
  
if(argc!=3){printf(" Use: ./corrections_calculator.exe hist_div_runXXXX.txt fit_t0_runXXXX.txt \n"); return 0;}

 TString namefile1= argv[1];
 TString namefile2= argv[2];

  int wheel[250],station[250],sector[250];
  float mean[250],meanerr[250],sigma[250],sigmaerr[250];
  float meant0[250],meanerrt0[250],sigmat0[250],sigmaerrt0[250],dummy[250];
  float correction[250],correctionerr[250];
  TH1F *correction_hist[4];
  TH1F *correction_all = new TH1F("correction","<t0>-worst_phase+12.5",32,-3,5);
  THStack *hs = new THStack("hs","");
  //TH1F *h_fits[6];
  //for(int a=0;a<6;a++) h_fits[a]= new TH1F(Form("peaks_histogram_region%d",a),Form("peaks_histogram_region%d",a),32,-3,5);

  for(int station=0;station<4;station++){
    correction_hist[station]= new TH1F(Form("Correction_for_MB%d",station+1),Form("Correction_for_MB%d",station+1),40,0,20);
  }
  
  ifstream in;
  in.open(namefile1);

  for(int row=0;row<250;row++ ){
    in>>wheel[row]>>station[row]>>sector[row]>>mean[row]>>meanerr[row]>>sigma[row]>>sigmaerr[row];
  }
  in.close();

  ifstream in2;
  in2.open(namefile2);

  for(int row=0;row<250;row++ ){
    in2>>wheel[row]>>station[row]>>sector[row]>>meant0[row]>>meanerrt0[row]>>sigmat0[row]>>sigmaerrt0[row]>>dummy[row];
  }
  in2.close();

  //////////////////////////////CALCULATE CORRECTIONS AND FILL/////////////////////////////////////////

  TString namefile3=namefile1.ReplaceAll("hist_div_","");
  TString Output=namefile3.ReplaceAll(".txt","");

  ofstream correctionf(Form("corrections_%s.txt",Output.Data()));
  correctionf<< " ...........Correction: <t0>-worst_phase_peak +12.5............................... "  <<endl;
  correctionf<<"    Wheel "<<"  "<<" MB "<<"  "<<" Sector "<<" "<< " correction "<<" correction error "<<endl;
  correctionf<< " ................................................................................. "  <<endl;
  int row=0;
  cout<<"       Generating corrections file "<<"correction_hist_"<<Output.Data()<<".txt ....."<<endl;
  for(int w=0;w<5;w++ ){
    for(int st=0;st<4;st++ ){
      for(int sec=0;sec<14;sec++ ){
	if(st!=3&&sec>11)continue;
	if(sec!=0&&sec!=6)correction[row]=meant0[row]-mean[row]+12.5;
	if(sec==6)correction[row]=meant0[row]-(mean[row-1]+mean[row+1])/2+12.5; ///To solve problems with vertical sectors in the cosmics analysis
	if(sec==0)correction[row]=meant0[row]-(mean[row+11]+mean[row+1])/2+12.5;
	correctionerr[row]=TMath::Sqrt(meanerrt0[row]*meanerrt0[row]+meanerr[row]*meanerr[row]);
	//cout<<"wheel "<<w-2<<" "<<" station "<<st+1<<" sector "<<sec+1<<endl;
	//cout<<mean[row]<<"  "<<meant0[row]<<"  "<<correction[row]<<endl;
	correction_hist[st]->Fill(correction[row]);
	correction_all->Fill(correction[row]);
	//Stacked coloured plot divided in 6 parts (June 2017)
	/*
	if(st==3&&sec<5&&sec>1)h_fits[4]->Fill(correction[row]);
	else{
	  if(st==3&&sec==12)h_fits[4]->Fill(correction[row]);
	  else{
	    if(st==0&&w!=1&&w!=2&&w!=3)h_fits[5]->Fill(correction[row]);
	    else h_fits[st]->Fill(correction[row]);
	  }
	}
	*/
	correctionf<<"      "<<showpos<<w-2<<noshowpos<<"      "<<st+1<<"      "<<setfill('0') << setw(2)<< sec+1<<"      "<< fixed << setprecision(4)<<correction[row]<<setw(2)<<"   "<<correctionerr[row]<<endl;	
	row++;
      }
    }
  }
  cout<<"  ........ Done!"<<endl;
  /////////////////////////////////////////CORRECTION PLOTS///////////////////////////////////////////
  
  TCanvas *canhist_corrections = new TCanvas(); 
  canhist_corrections  ->Divide(4,1);
  TF1 *myfitnew = new TF1("myfitnew","gaus(0)",0,15);
 
  gStyle->SetOptFit(0111);
  gStyle->SetOptStat(0);

  canhist_corrections->cd(1);
  /*
  for(int station=0;station<4;station++){
    canhist_corrections->cd(station+1);
    correction_hist[station] ->SetName(Form("MB%d",station+1));     
    correction_hist[station] ->Draw("");
    correction_hist[station] ->SetTitle(Form("MB%d",station+1));
    correction_hist[station] ->GetXaxis()->SetTitle("Correction [ns]");
    correction_hist[station] ->GetYaxis()->SetTitle("Entries");
    correction_hist[station] ->GetYaxis()->SetRangeUser(0,30);
    correction_hist[station] ->SetLineColor(station+1);
    int binmax  = correction_hist[station]->GetMaximumBin();
    double xmax = correction_hist[station]->GetXaxis()->GetBinCenter(binmax);
    double ymax = correction_hist[station]->GetMaximum();
    myfitnew->SetParameter(0,ymax);
    myfitnew->SetParameter(1,xmax);
    myfitnew->SetParameter(2,2);
    myfitnew->SetParName(0,"peak_height");
    myfitnew->SetParName(1,"peak_center");
    myfitnew->SetParName(2,"peak_sigma");
    correction_hist[station]->Fit("myfitnew","R");
  }
  canhist_corrections ->Print(Form("correction_hist_%s.gif",Output.Data()));
  */
  TCanvas *corrcanv= new TCanvas();
  correction_all->Draw();
  correction_all->SetTitle("correction: <t0>-worst_phase+12.5");
  correction_all->GetXaxis()->SetTitle("time(ns)");
  correction_all->GetYaxis()->SetTitle("entries");
  corrcanv->Print(Form("correction_all_%s.gif",Output.Data())); 
  corrcanv->Close();
  
 TCanvas *canstack= new TCanvas();
  for(int a=0;a<6;a++){
    h_fits[a]->SetFillColor(a+1);
    hs->Add(h_fits[a]);
  }  
  hs->Draw();
  hs->GetXaxis()->SetTitle("correction [ns]");
  hs->GetYaxis()->SetTitle("number of chambers");
  hs->SetTitle("corrections for the different regions");
  //canstack->Print("corrections_6regions.gif");
  canstack->Print("corrections.gif");


}
