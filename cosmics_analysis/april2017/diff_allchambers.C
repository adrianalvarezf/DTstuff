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
#include <string>
#include <vector>
#include <stdlib.h>

#endif

int main(int argc, const char* argv[])
{

  if(argc!=3){printf("Not enough arguments provided. Use: ./diff.exe file1.txt file2.txt\n"); return 0;}
  TString namefile1 =argv[1];
  TString namefile2 =argv[2];

  int wheel[250],station[250],sector[250];
  float mean[250],meanerr[250], sigma[250],sigmaerr[250];
  float mean2[250],meanerr2[250], sigma2[250],sigmaerr2[250];
  double difference[250],differenceerr[250];
  TH1F *difference_hist = new TH1F("difference","phase min",50,-1.5,1.5);

  std::ifstream in;
  in.open(namefile1);
  for(int row=0;row<250;row++ ){
    in>>wheel[row]>>station[row]>>sector[row]>>mean[row]>>meanerr[row]>>sigma[row]>>sigmaerr[row];
       std::cout<<wheel[row]<<"  "<<station[row]<<"  "<<sector[row]<<"  "<<mean[row]<<"  "<<meanerr[row]<<"  "<<sigma[row]<<"  "<<sigmaerr[row]<<std::endl;
  }
  in.close();

  std::ifstream in2;
  in2.open(namefile2);
  for(int row=0;row<250;row++ ){
    in2>>wheel[row]>>station[row]>>sector[row]>>mean2[row]>>meanerr2[row]>>sigma2[row]>>sigmaerr2[row];
    std::cout<<wheel[row]<<"  "<<station[row]<<"  "<<sector[row]<<"  "<<mean2[row]<<"  "<<meanerr2[row]<<"  "<<sigma2[row]<<"  "<<sigmaerr2[row]<<std::endl;
  }
  in2.close();
 
  TString file1=namefile1.ReplaceAll("hist_div_","");
  TString file2=namefile2.ReplaceAll("hist_div_","");
  TString runa=file1.ReplaceAll(".txt","");
  TString runb=file2.ReplaceAll(".txt","");

  std::ofstream txtdiff(Form("difference_%s_%s.txt",runa.Data(),runb.Data()));
  txtdiff<< " ........Difference: phase min "<<namefile1<<" and "<<namefile2<<"..." <<std::endl;
  txtdiff<< "    Wheel "<<"  "<<" MB "<<"  "<<" Sector "<<" "<< " difference "<<"  "<< " difference error "<<std::endl;
  txtdiff<< " ........................................................................... " <<std::endl;
  int row=0;

  for(int w=0;w<5;w++ ){
    for(int st=0;st<4;st++ ){
      for(int sec=0;sec<14;sec++ ){
	//if(st==1||st==2)continue;       //changed
	//if(w==1||w==3)continue;         //chambers
	if(sec>11&st!=3)continue;
	difference[row]=mean[row]-mean2[row];
	differenceerr[row]=std::sqrt(meanerr[row]*meanerr[row]+meanerr2[row]*meanerr2[row]);
	std::cout<<mean[row]<<"  "<<mean2[row]<<"  "<<" difference "<<difference[row]<<" difference error "<<differenceerr[row]<<std::endl;
	difference_hist->Fill(difference[row]);
	txtdiff<<"      "<<std::showpos<<w-2<<std::noshowpos<<"      "<<st+1<<"      "<<std::setfill('0') << std::setw(2)<< sec+1<<"      "<< std::fixed <<std::showpos<< difference[row]<<std::noshowpos<<"    "<<differenceerr[row]<<std::endl;	
	row++;
      }
    }
  }

  TCanvas *diffcanv= new TCanvas();
  difference_hist->Draw();
  difference_hist->SetTitle(Form("Difference in phase min %s and %s",runa.Data(),runb.Data()));
  difference_hist->GetXaxis()->SetTitle("time(ns)");
  difference_hist->GetYaxis()->SetTitle("entries");
  diffcanv->Print(Form("difference_%s_%s.gif",runa.Data(),runb.Data())); 
    
  return 0;
 
}
