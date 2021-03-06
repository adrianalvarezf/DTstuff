
#define t0_fitter_chain_cxx

#include "t0_fitter_chain.h"
#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <string>
#include "Riostream.h"
#include "TVectorT.h"
#include "TString.h"
#include <iostream>

///////////////// Version 6 has 14 sectors//////////////////////////


void t0_fitter_chain::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

Bool_t t0_fitter_chain::Process(Long64_t entry)
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
  vector<short>   &my_ltTwinMuxIn_bx= *ltTwinMuxIn_bx;
 
  vector<float>   &my_dtsegm4D_x_dir_loc = *dtsegm4D_x_dir_loc;
  vector<float>   &my_dtsegm4D_z_dir_loc = *dtsegm4D_z_dir_loc;

  vector<short>   &my_Mu_isMuGlobal = *Mu_isMuGlobal;
  vector<short>   &my_Mu_isMuTracker = *Mu_isMuTracker;
  vector<float>   &my_Mu_normchi2_glb = *Mu_normchi2_glb;
  vector<int>     &my_Mu_numberOfMatches_sta = *Mu_numberOfMatches_sta;
  vector<float>   &my_Mu_dz_glb = *Mu_dz_glb;
  vector<float>   &my_Mu_dxy_glb = *Mu_dxy_glb;
  vector<int>     &my_Mu_numberOfPixelHits_glb = *Mu_numberOfPixelHits_glb;
  vector<int>     &my_Mu_numberOfTrackerHits_glb = *Mu_numberOfTrackerHits_glb;
  vector<int>     &my_Mu_numberOfHits_sta = *Mu_numberOfHits_sta;
  vector<int>     &my_Mu_segmentIndex_sta = *Mu_segmentIndex_sta;
  vector <float>  &my_Mu_px = *Mu_px;
  vector <float>  &my_Mu_py = *Mu_py;
  vector <float>  &my_Mu_pz = *Mu_pz;
  vector <float>  &my_Mu_phi = *Mu_phi;
  vector <float>  &my_Mu_eta = *Mu_eta;

 
  //fChain-> Add ("/eos/cms/store/group/dpg_dt/comm_dt/PromptOffline/POATrees2017/ExpressPhysics/crab_runTree_303832/170925_124218/0000/DTNtuple_2.root",-1);
  //fChain-> Add ("/eos/cms/store/group/dpg_dt/comm_dt/PromptOffline/POATrees2017/ExpressPhysics/crab_runTree_303832/170925_124218/0000/DTNtuple_3.root",-1);
  fChain->GetTree()->GetEntry(entry);
 
  if(fmod(entry,10000)==0) Printf(" ..... event %d", int(entry));


  int ndtsegm4D =  dtsegm4D_station->size();
  
  /*
    nseghist->Fill(ndtsegm4D);
    if(Nmuons ==1) nseghist_1muon->Fill(ndtsegm4D);
    if(Nmuons ==1 && my_Mu_isMuGlobal[0]==1 && my_Mu_isMuTracker[0]==1 && my_Mu_normchi2_glb[0]<10&& my_Mu_numberOfHits_sta[0]>0 && my_Mu_numberOfMatches_sta[0]>1 && fabs(my_Mu_dxy_glb[0])<0.2 && my_Mu_numberOfPixelHits_glb[0]>0 && my_Mu_numberOfTrackerHits_glb[0]>5){
    nseghist_wcuts->Fill(ndtsegm4D);}
  */
  vector<int>  used_segments;
  used_segments.clear();
  /*
  for (int idtsegm = 0; idtsegm < ndtsegm4D; ++idtsegm) {
    float ang_rec=TMath::ATan(my_dtsegm4D_x_dir_loc[idtsegm]/my_dtsegm4D_z_dir_loc[idtsegm]);
    TVectorF *my_segm4D_phi_hitsLayer = (TVectorF*)dtsegm4D_phi_hitsLayer->At(idtsegm);
    int n_phi_hitsLayer = my_segm4D_phi_hitsLayer->GetNoElements();
    if(n_phi_hitsLayer<6)continue;
    //if(fabs(ang_rec) >TMath::Pi()/6)continue;
    ang_hist->Fill(ang_rec*180/(TMath::Pi()));
  }
  */
  
  mtight=0;
  for(int nmu=0;nmu<Nmuons;nmu++){
    if(my_Mu_isMuGlobal[nmu]!=1)continue;
    if(my_Mu_isMuTracker[nmu]!=1)continue;
    if(my_Mu_normchi2_glb[nmu]>10)continue;
    if(my_Mu_numberOfHits_sta[nmu]<1)continue;
    if(my_Mu_numberOfMatches_sta[nmu]<2)continue;
    if(fabs(my_Mu_dxy_glb[nmu])>0.2)continue;
    if(my_Mu_numberOfPixelHits_glb[nmu]<=0)continue;
    if(my_Mu_numberOfTrackerHits_glb[nmu]<=5)continue;
    mtight++;
    int segmWord = my_Mu_segmentIndex_sta[nmu];
    for (int idtsegm = 0; idtsegm < ndtsegm4D; ++idtsegm) {
      if(my_ltTwinMuxIn_quality.size()==0)continue;
       int bit_seg = (segmWord >> idtsegm) & 1;
      //cout << " Segment " << idtsegm << " bit " << bit_seg << endl; 
      TVectorF *my_segm4D_phi_hitsLayer = (TVectorF*)dtsegm4D_phi_hitsLayer->At(idtsegm);
      int n_phi_hitsLayer = my_segm4D_phi_hitsLayer->GetNoElements();
      float ang_rec=TMath::ATan(my_dtsegm4D_x_dir_loc[idtsegm]/my_dtsegm4D_z_dir_loc[idtsegm]);
      if(fabs(ang_rec) >TMath::Pi()/6)continue;
      if(n_phi_hitsLayer<6)continue;
      // if(TMath::Sqrt(my_Mu_px[nmu]*my_Mu_px[nmu]+my_Mu_py[nmu]*my_Mu_py[nmu])<25)continue;
      if(bit_seg!=1)continue;
       int used=0;
       for(int seg=0;seg<used_segments.size();seg++){
	 segmentstotal++;
	 if(used_segments[seg]==idtsegm){
	   used=1; 
	   segmentsfortwo++;
	 }
       }
       if(used>0)continue;
     
       for (int q = 0; q <(int)my_ltTwinMuxIn_quality.size(); ++q) {
	if(my_ltTwinMuxIn_station[q]!= my_dtsegm4D_station[idtsegm] || my_ltTwinMuxIn_wheel[q]!= my_dtsegm4D_wheel[idtsegm] ){continue;}
	if(my_ltTwinMuxIn_sector[q]!= my_dtsegm4D_sector[idtsegm] &&  my_dtsegm4D_sector[idtsegm]<13 ){continue;}
	if(my_dtsegm4D_sector[idtsegm]==13 && my_ltTwinMuxIn_sector[q]!=4){continue;}
	if(my_dtsegm4D_sector[idtsegm]==14 && my_ltTwinMuxIn_sector[q]!=10){continue;}
	if(fabs(my_ltTwinMuxIn_bx[q])>1)continue;
	T0[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill(my_dtsegm4D_t0[idtsegm]);
	hist_bx->Fill(my_ltTwinMuxIn_bx[q]);
	mu_phi->Fill(my_Mu_phi[nmu]);
	mu_eta->Fill(my_Mu_eta[nmu]);
	ang_hist->Fill(ang_rec*180/(TMath::Pi()));
     	mu_pt->Fill(TMath::Sqrt(my_Mu_px[nmu]*my_Mu_px[nmu]+my_Mu_py[nmu]*my_Mu_py[nmu]));
	BX[my_ltTwinMuxIn_wheel[q]+2][my_ltTwinMuxIn_station[q]-1][my_dtsegm4D_sector[idtsegm]-1]->Fill(my_ltTwinMuxIn_bx[q]);
      }
      used_segments.push_back(idtsegm);
    }
  }
  mu_tight->Fill(mtight);
 
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

  my_Mu_isMuGlobal.clear();
  my_Mu_isMuTracker.clear();
  my_Mu_normchi2_glb.clear();
  my_Mu_numberOfMatches_sta.clear();
  my_Mu_dxy_glb.clear();
  my_Mu_dz_glb.clear();
  my_Mu_numberOfHits_sta.clear();
  my_Mu_numberOfPixelHits_glb.clear();
  my_Mu_numberOfTrackerHits_glb.clear();
  my_Mu_px.clear();
  my_Mu_py.clear();
  my_Mu_pz.clear(); 
  my_Mu_phi.clear();
  my_Mu_eta.clear(); 

  return kTRUE;

  }
 
void t0_fitter_chain::Terminate()
{

 
  TFile *my_new_file = new TFile(Form("t0_histograms_run%d.root",runnumber),"RECREATE"); 
  //TFile *my_new_file = new TFile("t0_histograms_Run295463_and_296173.root","RECREATE"); 
  TCanvas *can3[5];
  TGraphErrors *gr[5][4];
  int sec[14]={1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      gr[wheel][station]= new TGraphErrors();
      for(int sector=0;sector<14;sector++){
	double mean= T0[wheel][station][sector]->GetMean();
	double sigma= T0[wheel][station][sector]->GetRMS();
	gr[wheel][station]->SetPoint(sector,sec[sector]+0.2*(station+1),mean);
	gr[wheel][station]->SetPointError(sector,0,sigma);
      } 
    }
  }
 
  for(int wheel=0;wheel<5;wheel++){
    can3[wheel]= new TCanvas(); 
    can3[wheel]->cd(1);
    TString namet0 = Form("t0_Wh%d_average.gif",wheel-2);  
    for(int station=0;station<4;station++){
      gr[wheel][station]->SetName(Form("t0_Wh%d_st%d.gif",wheel-2,station+1));     
      if(station==0){
	gr[wheel][station]->Draw("AP");
	gr[wheel][station]->SetTitle(namet0);
	gr[wheel][station]->GetXaxis()->SetTitle("Sector number");
	gr[wheel][station]->GetYaxis()->SetTitle("<t0>");
	gr[wheel][station]->GetYaxis()->SetRangeUser(-10.,10.);
	gr[wheel][station]->GetXaxis()->SetRangeUser(0.,15.);
      }
      else gr[wheel][station]->Draw("SAMEP");
      gr[wheel][station]->SetMarkerColor(station+1);
      gr[wheel][station]->SetMarkerStyle(station+20);
      gr[wheel][station]->SetMarkerSize(1.);
      gr[wheel][station]->SetLineColor(station+1);
    }
    can3[wheel]->Print(namet0);  
    can3[wheel]->Write();  
    can3[wheel]->Close();
  }
  

  TCanvas *can4[5][4][14];
  //ofstream fitdata ("fit_t0_Run295463_and_296173.txt");
  ofstream fitdata (Form("fit_t0_run%d.txt",runnumber));
 
  fitdata<<"    Wheel "<<"  "<<" MB "<<" "<<" Sector "<<"    "<< " mean "<<" "<<" mean error "<<" "<<" sigma "<<" sigma error "<<"   "<<" BX "<<endl;
  fitdata<< " ................................................................................. " <<endl;

  float mainpeak=0;
  float secondpeak=0;
  float thirdpeak=0;
  
  for(int w=0;w<5;w++){
    for(int st=0;st<4;st++){
      for(int sec=0;sec<14;sec++){
	if(sec>11 && st!=3)continue;
	if(T0[w][st][sec]->GetEntries()==0){
	  fitdata<<"missing chamber Wh"<<w-2<<" MB"<<st+1<<" Sec"<<sec+1<<endl; 
	  continue;
	}
	TString namet0s = Form("t0_Wh%d_MB%d_Sec%d.gif",w-2,st+1,sec+1); 
  	can4[w][st][sec]= new TCanvas(); 
	can4[w][st][sec]->cd(1);
 	T0[w][st][sec]->Draw();
	T0[w][st][sec]->Fit("gaus","","",-10,10);
	T0[w][st][sec]->SetTitle(namet0s);
	T0[w][st][sec]->GetXaxis()->SetTitle("t0");
	T0[w][st][sec]->GetYaxis()->SetTitle("entries");
	TF1 *fit1 = (TF1*)T0[w][st][sec]->GetFunction("gaus");
	double mean = fit1->GetParameter(1);
	double sigma = fit1->GetParameter(2);
	double meanerr = fit1->GetParError(1);
	double sigmaerr = fit1->GetParError(2);
	myfit->SetParName(0,"height");
	myfit->SetParName(1,"mean");
	myfit->SetParName(2,"sigma");
	gStyle->SetOptFit(0111);
	gStyle->SetOptStat(0);
	double bxmean= BX[w][st][sec]->GetMean();
	h_fits[st]->Fill(mean);
	if(st==3 && (sec==2||sec==3||sec==4||sec==12))thirdpeak+=mean;
	else if (st==0 && (w==0 ||w==4))secondpeak+=mean;
	else mainpeak+=mean;
	fitdata<<"      "<<showpos<<w-2<<noshowpos<<"      "<<st+1<<"      "<<setfill('0') << setw(2)<< sec+1<<"      "<< fixed << setprecision(4)<<showpos<<mean<<"   "<<noshowpos<<meanerr<<"     "<<fabs(sigma)<<"   "<<sigmaerr<<"     "<<showpos<<bxmean<<endl;
	//printf("wheel %3d station %2d sector %3d <t0>= %7.4f +/- %7.4f sigma= %7.4f +/- %7.4f\n",w-2,st+1,sec+1,mean,meanerr,sigma,sigmaerr);
	can4[w][st][sec]->Print(namet0s); 
	T0[w][st][sec]->Write();
	can4[w][st][sec]->Close(); 
      }
    }
  }
  secondpeak=secondpeak/28;
  thirdpeak=thirdpeak/20;
  mainpeak=mainpeak/202;
  cout<<"main peak=   "<<mainpeak<<"   second peak=   "<<secondpeak<<"    third peak=   "<<thirdpeak<<endl;
  TCanvas *histomean= new TCanvas();
  for(int a=0;a<4;a++){
    h_fits[a]->SetFillColor(a+1);
    hs->Add(h_fits[a]);
  }  
  hs->Draw();
  hs->GetXaxis()->SetTitle("<t0> [ns]");
  hs->GetYaxis()->SetTitle("entries");
  hs->SetTitle("<t0> from fits");
  histomean->Print("histo_t0means.gif");

  
  /*
  TFile *my_new_file2 = new TFile("histograms_of_muon_variables.root","RECREATE"); 
  TCanvas *angcanv =new TCanvas();
  ang_hist->Draw();
  ang_hist->GetYaxis()->SetTitle("entries");
  ang_hist->GetXaxis()->SetTitle("segment angle (degrees)");
  angcanv->Print("angseg_hist.gif");
  angcanv->Close();
  ang_hist->Write();
  TCanvas *phicanv =new TCanvas();
  mu_phi->Draw();
  mu_phi->GetYaxis()->SetTitle("entries");
  mu_phi->GetXaxis()->SetTitle("muon angle phi in rad");
  phicanv->Print("muonphi_hist.gif");
  phicanv->Close();
  mu_phi->Write();
  TCanvas *etacanv =new TCanvas();
  mu_eta->Draw();
  mu_eta->GetYaxis()->SetTitle("entries");
  mu_eta->GetXaxis()->SetTitle("muon pseudorapidity");
  etacanv->Print("muoneta_hist.gif");
  etacanv->Close();
  mu_eta->Write();
  TCanvas *ptcanv =new TCanvas();
  mu_pt->Draw();
  mu_pt->GetYaxis()->SetTitle("entries");
  mu_pt->GetXaxis()->SetTitle("muon pt");
  ptcanv->Print("muonpt_hist.gif");
  ptcanv->Close();
  mu_pt->Write();
  TCanvas *bxcanv =new TCanvas();
  hist_bx->Draw();
  hist_bx->GetYaxis()->SetTitle("entries");
  hist_bx->GetXaxis()->SetTitle("BX");
  bxcanv->Print("bx_hist.gif");
  bxcanv->Close();
  hist_bx->Write();
  TCanvas *tightcanv =new TCanvas();
  mu_tight->Draw();
  mu_tight->GetYaxis()->SetTitle("entries");
  mu_tight->GetXaxis()->SetTitle("Number of tight muons");
  tightcanv->Print("tightmu.gif");
  tightcanv->Close();
  tightcanv->Write();
  */

  //cout<<" segments for two muons = "<<segmentsfortwo<<" segments total = "<<segmentstotal<<" fraction "<<(double)segmentsfortwo/segmentstotal<<endl;

  
}
