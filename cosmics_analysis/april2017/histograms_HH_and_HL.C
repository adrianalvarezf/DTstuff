
#define histograms_HH_and_HL_cxx

#include "histograms_HH_and_HL.h"
#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <string>
#include "Riostream.h"
#include "TVectorT.h"
#include "TString.h"

void histograms_HH_and_HL::Begin(TTree * /*tree*/)
{
 
  TString option = GetOption();

}

Bool_t histograms_HH_and_HL::Process(Long64_t entry)
{
  vector<short>   &my_dtsegm4D_wheel= *dtsegm4D_wheel;
  vector<short>   &my_dtsegm4D_station= *dtsegm4D_station; 
  vector<short>   &my_dtsegm4D_sector= *dtsegm4D_sector;
  vector<short>   &my_dtsegm4D_phinhits= *dtsegm4D_phinhits;
  vector<float>   &my_dtsegm4D_t0 = *dtsegm4D_t0;
  vector<short>   &my_ltTwinMuxIn_quality= *ltTwinMuxIn_quality;
  vector<short>   &my_ltTwinMuxIn_wheel= *ltTwinMuxIn_wheel;
  vector<short>   &my_ltTwinMuxIn_sector= *ltTwinMuxIn_sector;
  vector<short>   &my_ltTwinMuxIn_station= *ltTwinMuxIn_station;
  vector<float>   &my_dtsegm4D_x_dir_loc = *dtsegm4D_x_dir_loc;
  vector<float>   &my_dtsegm4D_z_dir_loc = *dtsegm4D_z_dir_loc;
s
  Long64_t nentries = fChain->GetEntriesFast();
  fChain->GetEntry(entry);
  if(fmod(entry,10000)==0) Printf(" ..... event %d of %d", int(entry), int(nentries));


  int ndtsegm4D =  dtsegm4D_station->size();
  totalsegm+=ndtsegm4D;
  ////////////////////////////////////////////////////////////////////FILLING HISTOS////////////////////////////////////////////////////////////////
  for (int idtsegm = 0; idtsegm < ndtsegm4D; ++idtsegm) {
    if(my_ltTwinMuxIn_quality.size()==0){continue;}
    TVectorF *my_segm4D_phi_hitsLayer = (TVectorF*)dtsegm4D_phi_hitsLayer->At(idtsegm);
    int n_phi_hitsLayer = my_segm4D_phi_hitsLayer->GetNoElements();
    float ang_rec=TMath::ATan(my_dtsegm4D_x_dir_loc[idtsegm]/my_dtsegm4D_z_dir_loc[idtsegm]);
    if(fabs(ang_rec) >TMath::Pi()/6)continue;
    //if(n_phi_hitsLayer<6){continue;}

    ///TwinMux and dtsegm4D chambers should match
    for (int q = 0; q <(int)my_ltTwinMuxIn_quality.size(); ++q) {
      if(my_ltTwinMuxIn_station[q]!= my_dtsegm4D_station[idtsegm] || my_ltTwinMuxIn_wheel[q]!= my_dtsegm4D_wheel[idtsegm]){continue;}
      if(my_ltTwinMuxIn_sector[q]!= my_dtsegm4D_sector[idtsegm] &&  my_dtsegm4D_sector[idtsegm]<13){continue;}
      if(my_dtsegm4D_sector[idtsegm]==13 && my_ltTwinMuxIn_sector[q]!=4){continue;}
      if(my_dtsegm4D_sector[idtsegm]==14 && my_ltTwinMuxIn_sector[q]!=10){continue;}
      if(my_ltTwinMuxIn_quality[q]==6){
	HH[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill(my_dtsegm4D_t0[idtsegm]);
	HH_fold[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill((int)(my_dtsegm4D_t0[idtsegm]+50)%25);
	highhigh++;both++;
	break;
      }
      if(my_ltTwinMuxIn_quality[q]==5){
	HL[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill(my_dtsegm4D_t0[idtsegm]);
	HL_fold[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill((int)(my_dtsegm4D_t0[idtsegm]+50)%25);
	highlow++;both++;
	break;
      }
      if(my_ltTwinMuxIn_quality[q]<5){worseq++;	break;}
    }
  }
 
 
 return kTRUE;

 my_dtsegm4D_wheel.clear();
 my_dtsegm4D_station.clear(); 
 my_dtsegm4D_sector.clear();
 my_dtsegm4D_phinhits.clear();
 my_dtsegm4D_t0.clear();
 my_ltTwinMuxIn_quality.clear();
 my_ltTwinMuxIn_sector.clear();
 my_ltTwinMuxIn_station.clear();
 my_ltTwinMuxIn_wheel.clear();
 my_dtsegm4D_x_dir_loc.clear();
 my_dtsegm4D_z_dir_loc.clear();


}
 
void histograms_HH_and_HL::Terminate()
{

  cout<<" Number of HH = "<<highhigh<<"," <<" Number of HL = "<<highlow<<"," <<" Number of both = "<<both<<","<<" Number of bad quality segments "<<worseq<<","<<" Number of segments "<<totalsegm<<endl;
  
  //TFile *my_new_file = new TFile(Form("run%d_%dV_histograms.root",runnumber,voltage),"RECREATE"); 
  TFile *my_new_file = new TFile(Form("run%d_histograms.root",runnumber),"RECREATE"); 

  ////////////////////////////////////////////////////////////////////NON-FOLDED HISTOGRAMS////////////////////////////////////////////////////////////////
  TCanvas *can2[5][4][14][2];
  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sec=0;sec<14;sec++){
	for(int qu=0;qu<2;qu++){
	  TString name[2]={Form("HL_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1),Form("HH_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1)};
	  if(sec>11&&station!=3)continue;
	  can2[wheel][station][sec][qu]= new TCanvas(); 
	  can2[wheel][station][sec][qu]->cd(1);
	  if(qu==1){
	    HH[wheel][station][sec]->Draw();
	    HH[wheel][station][sec]->Draw("E1SAME");
	    HH[wheel][station][sec]->SetBins(120,-80,40);
	    HH[wheel][station][sec]->SetTitle(name[qu]);
	    HH[wheel][station][sec]->GetXaxis()->SetTitle("t0(ns)");
	    HH[wheel][station][sec]->GetYaxis()->SetTitle("entries");
	    HH[wheel][station][sec]->Write();
	  }
	  if(qu==0){
	    HL[wheel][station][sec]->Draw();
	    HL[wheel][station][sec]->Draw("E1SAME");
	    HL[wheel][station][sec]->SetBins(120,-80,40);
	    HL[wheel][station][sec]->SetTitle(name[qu]);
	    HL[wheel][station][sec]->GetXaxis()->SetTitle("t0(ns)");
	    HL[wheel][station][sec]->GetYaxis()->SetTitle("entries");
	    HL[wheel][station][sec]->Write();
	  }
	  //can2[wheel][station][sec][qu]->Print(name[qu]);  //Uncomment to see non-folded histograms
	  can2[wheel][station][sec][qu]->Close();  
	} 
      }
    }
  }

////////////////////////////////////////////////////////////////////FOLDED HISTOGRAMS////////////////////////////////////////////////////////////////
 TCanvas *can22[5][4][14][2];
  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sec=0;sec<14;sec++){
	for(int qu=0;qu<2;qu++){
	  TString name[2]={Form("HL_fold_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1),Form("HH_fold_Wh%d_MB%d_Sec%d.gif",wheel-2,station+1,sec+1)};
	  if(sec>11&&station!=3)continue;
	  can22[wheel][station][sec][qu]= new TCanvas(); 
	  can22[wheel][station][sec][qu]->cd(1);
	  if(qu==1){
	    HH_fold[wheel][station][sec]->Draw();
	    HH_fold[wheel][station][sec]->Draw("E1SAME");
	    HH_fold[wheel][station][sec]->SetBins(50,0,25);
	    HH_fold[wheel][station][sec]->SetTitle(name[qu]);
	    HH_fold[wheel][station][sec]->GetXaxis()->SetTitle("t0(ns)");
	    HH_fold[wheel][station][sec]->GetYaxis()->SetTitle("entries");
	    HH_fold[wheel][station][sec]->Write();	 
	  }
	  if(qu==0){
	    HL_fold[wheel][station][sec]->Draw();
	    HL_fold[wheel][station][sec]->Draw("E1SAME");
	    HL_fold[wheel][station][sec]->SetBins(50,0,25);
	    HL_fold[wheel][station][sec]->SetTitle(name[qu]);
	    HL_fold[wheel][station][sec]->GetXaxis()->SetTitle("t0(ns)");
	    HL_fold[wheel][station][sec]->GetYaxis()->SetTitle("entries");
	    HL_fold[wheel][station][sec]->Write();
	  }
	  can22[wheel][station][sec][qu]->Print(name[qu]);  
	  can22[wheel][station][sec][qu]->Close();  
	} 
      }
    }
  }
  

  

  
}
