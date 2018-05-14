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

#endif


int main(int argc, const char* argv[]){

  if(argc=!4){printf("Not enough arguments provided. Use: ./histodiv.exe [numberoffiles N] [print gifs 1=y/0=n] file1.root ...  fileN.root\n"); return 0;}

  //////////////////////////////////////////////////////////DEFINITIONS//////////////////////////////////////////////////
  TH1F *HH[5][4][14];  
  TH1F *HL[5][4][14]; 
  TH1F *HH_clone[5][4][14];  
  TH1F *HL_clone[5][4][14]; 
  TH1F *div[5][4][14];
  TH1F *HH_fold[5][4][14];  
  TH1F *HL_fold[5][4][14];
  TH1F *HH_clone_fold[5][4][14];  
  TH1F *HL_clone_fold[5][4][14];
  TH1F *div_fold[5][4][14]; 
  TF1  *myfit = new TF1("myfit","gaus(0)+pol0(3)",5,15);
  int filenum= atoi(argv[1]);
  int printgif= atoi(argv[2]);

  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sector=0;sector<14;sector++){
	TString hName = Form("Wh%d_St%d_Sec%d",wheel-2,station+1,sector+1);
	TString hNamef = Form("Wh%d_St%d_Sec%d",wheel-2,station+1,sector+1);
	
	HL[wheel][station][sector] = new TH1F(Form("%s_HL",hName.Data()),hName,120,-80,40); 
	HH[wheel][station][sector] = new TH1F(Form("%s_HH",hName.Data()),hName,120,-80,40); 
	HL_fold[wheel][station][sector] = new TH1F(Form("%s_HL_fold",hNamef.Data()),hName,50,0,25); 
	HH_fold[wheel][station][sector] = new TH1F(Form("%s_HH_fold",hNamef.Data()),hName,50,0,25); 
	 
	HL_clone[wheel][station][sector] = new TH1F(Form("%s_HL_clone",hName.Data()),hName,120,-80,40); 
	HH_clone[wheel][station][sector] = new TH1F(Form("%s_HH_clone",hName.Data()),hName,120,-80,40); 
	HL_clone_fold[wheel][station][sector] = new TH1F(Form("%s_HL_clone_fold",hNamef.Data()),hName,50,0,25); 
	HH_clone_fold[wheel][station][sector] = new TH1F(Form("%s_HH_clone_fold",hNamef.Data()),hName,50,0,25); 

	div[wheel][station][sector]= new TH1F(Form("%s_div",hName.Data()),hName,120,-80,40);
	div_fold[wheel][station][sector]= new TH1F(Form("%s_div_fold",hNamef.Data()),hName,50,0,25); 
      }
    }
  }
  //////////////////////////////////////////////////////////HISTOGRAM FILLING//////////////////////////////////////////////////
  TFile *f[filenum];
  //f[0] = new TFile("run290910_3400V_histograms.root","READ");
  for(int i=0;i<filenum;i++) f[i] = new TFile(argv[i+3],"READ");

  for (int nfile=0;nfile<filenum;nfile++){
    for(int wheel=0;wheel<5;wheel++){
      for(int station=0;station<4;station++){
	for(int sec=0;sec<14;sec++){
	  if(sec>11&&station!=3)continue;
	  HH[wheel][station][sec] = (TH1F*)f[nfile]->Get(Form("Wh%d_St%d_Sec%d_HH",wheel-2,station+1,sec+1));
	  HL[wheel][station][sec] = (TH1F*)f[nfile]->Get(Form("Wh%d_St%d_Sec%d_HL",wheel-2,station+1,sec+1));
	  HH_fold[wheel][station][sec] = (TH1F*)f[nfile]->Get(Form("Wh%d_St%d_Sec%d_HH_fold",wheel-2,station+1,sec+1));
	  HL_fold[wheel][station][sec] = (TH1F*)f[nfile]->Get(Form("Wh%d_St%d_Sec%d_HL_fold",wheel-2,station+1,sec+1));
	  if(nfile==0){
	    HH_clone[wheel][station][sec]=(TH1F*)HH[wheel][station][sec]->Clone(); 
	    HL_clone[wheel][station][sec]=(TH1F*)HL[wheel][station][sec]->Clone(); 
	    HH_clone_fold[wheel][station][sec]=(TH1F*)HH_fold[wheel][station][sec]->Clone(); 
	    HL_clone_fold[wheel][station][sec]=(TH1F*)HL_fold[wheel][station][sec]->Clone(); 
	  }
	  else{
	    HH_clone[wheel][station][sec]->Add(HH[wheel][station][sec]);
	    HL_clone[wheel][station][sec]->Add(HL[wheel][station][sec]);
	    HH_clone_fold[wheel][station][sec]->Add(HH_fold[wheel][station][sec]);
	    HL_clone_fold[wheel][station][sec]->Add(HL_fold[wheel][station][sec]);
	  }
	}
      }
    }
  }
   
  TString sInput =argv[3];
  TString sBlabla=sInput.ReplaceAll("_histograms.root","");
  TString sOutput=Form("hist_div_%s",sBlabla.Data()); 
  TFile *my_new_file = new TFile(Form("%s.root",sOutput.Data()),"RECREATE");

  std::ofstream foldedfit (Form("%s.txt",sOutput.Data()));
  foldedfit<<"    Wheel "<<"  "<<" MB "<<" "<<" Sector "<<"    "<<" mean "<<" "<<" mean error "<<" "<<" sigma "<<" sigma error "<<std::endl;
  foldedfit<< " ................................................................................. " <<std::endl;

 
  TCanvas *can[5][4][14][2];
  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sec=0;sec<14;sec++){
	if(HH_clone[wheel][station][sec]->GetEntries()==0 || HL_clone[wheel][station][sec]->GetEntries()==0)continue;
	TString namediv = Form("HL_vs_HH_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1);
	TString namedivf = Form("HL_vs_HH_folded_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1);

	//////////////////////////////////////////////NON-FOLDED HISTOS///////////////////////////////////////////////////////////////////
	can[wheel][station][sec][0]= new TCanvas(); 
	can[wheel][station][sec][0]->cd(1);
	div[wheel][station][sec]=(TH1F*)HL_clone[wheel][station][sec]->Clone();
	div[wheel][station][sec]->Draw();    
	div[wheel][station][sec]->Draw("E1SAME");   
	for(int b=1;b<121;b++){
	  double hh= HH_clone[wheel][station][sec]->GetBinContent(b);
	  double hl= HL_clone[wheel][station][sec]->GetBinContent(b);
	  double err= TMath::Sqrt((TMath::Sqrt(hl)/hh)*(TMath::Sqrt(hl)/hh)+(TMath::Sqrt(hh)*hl/(hh*hh))*(TMath::Sqrt(hh)*hl/(hh*hh)));
	  div[wheel][station][sec]->SetBinError(b,err);
	} 
	div[wheel][station][sec]->Divide(HH_clone[wheel][station][sec]);
	div[wheel][station][sec]->GetYaxis()->SetRangeUser(0,6);
	div[wheel][station][sec]->GetXaxis()->SetTitle("t0(ns)");
	div[wheel][station][sec]->GetYaxis()->SetTitle("HL/HH");
	div[wheel][station][sec]->SetBins(120,-80,40);
	div[wheel][station][sec]->SetLineColor(2);
	div[wheel][station][sec]->SetTitle(namediv);
	//can[wheel][station][sec][0]->Print(namediv); // commented, no need to show non-folded histograms
 	gStyle->SetOptStat(0);
	can[wheel][station][sec][0]->Close();  

	//////////////////////////////////////////////FOLDED HISTOS + FITS + TXT OUTPUT///////////////////////////////////////////////////////////////////
	can[wheel][station][sec][1]= new TCanvas(); 
	can[wheel][station][sec][1]->cd(1);
	div_fold[wheel][station][sec]=(TH1F*)HL_clone_fold[wheel][station][sec]->Clone();
	div_fold[wheel][station][sec]->SetName(Form("HLvsHH_Wh%d_MB%d_Sec%d",wheel-2,station+1,sec+1));
	div_fold[wheel][station][sec]->Draw();    
	div_fold[wheel][station][sec]->Draw("E1SAME");   
	for(int bu=1;bu<51;bu++){
	  double hh= HH_clone_fold[wheel][station][sec]->GetBinContent(bu);
	  double hl= HL_clone_fold[wheel][station][sec]->GetBinContent(bu);
	  double err= TMath::Sqrt((TMath::Sqrt(hl)/hh)*(TMath::Sqrt(hl)/hh)+(TMath::Sqrt(hh)*hl/(hh*hh))*(TMath::Sqrt(hh)*hl/(hh*hh)));
	  div_fold[wheel][station][sec]->SetBinError(bu,err);
	}
	div_fold[wheel][station][sec]->Divide(HH_clone_fold[wheel][station][sec]);
	//div_fold[wheel][station][sec]->GetYaxis()->SetRangeUser(0,2.5);
	div_fold[wheel][station][sec]->GetXaxis()->SetTitle("t0_folded(ns)");
	div_fold[wheel][station][sec]->GetYaxis()->SetTitle("HL/HH");
	div_fold[wheel][station][sec]->SetBins(50,0,25);
	div_fold[wheel][station][sec]->SetLineColor(2);
	 
	int  binmax = div_fold[wheel][station][sec]->GetMaximumBin();
	double xmax = div_fold[wheel][station][sec]->GetXaxis()->GetBinCenter(binmax);
	double ymax = div_fold[wheel][station][sec]->GetMaximum();
	double ymin = div_fold[wheel][station][sec]->GetMinimum();
	myfit->SetParameter(0,ymax-ymin);
	myfit->SetParameter(1,xmax);
	myfit->SetParameter(2,3);
	//if(wheel==0&&station==0&&sec==8)myfit->SetParameter(1,5);
	//if(wheel==1&&station==0&&sec==0)myfit->SetParameter(1,9);
	myfit->SetParameter(3,ymin);
	div_fold[wheel][station][sec]->Fit("myfit","","",0,25);
	double finalmean = myfit->GetParameter(1);
	double finalsigma = myfit->GetParameter(2);
	double finalmeanerr = myfit->GetParError(1);
	double finalsigmaerr = myfit->GetParError(2);
	div_fold[wheel][station][sec]->SetTitle(namedivf);
	myfit->SetParName(0,"peak_height");
	myfit->SetParName(1,"peak_center");
	myfit->SetParName(2,"peak_sigma");
	myfit->SetParName(3,"initial_height");
	gStyle->SetOptFit(0111);
	gStyle->SetOptStat(0);
	TString namechamber = Form("Wh%d_MB%d_Sec%d",wheel-2,station+1,sec+1);
	foldedfit<<"      "<<std::showpos<<wheel-2<<std::noshowpos<<"      "<<station+1<<"      "<<std::setfill('0') << std::setw(2)<< sec+1<<"      "<< std::fixed << std::setprecision(4)<<std::setw(7)<<finalmean<<std::setw(2)<<"   "<<finalmeanerr<<"     "<<fabs(finalsigma)<<"   "<<finalsigmaerr<<std::endl;		
	div_fold[wheel][station][sec]->Write();
	if(printgif==1)can[wheel][station][sec][1]->Print(namedivf);  
	can[wheel][station][sec][1]->Close();  
	
      } 
    }
  }
  
  return 0;

}
