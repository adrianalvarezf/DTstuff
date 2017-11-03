#define MyEffWithDigis_all_cxx
#include "MyEffWithDigis_all.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TVectorF.h>
#include <TH2.h>
#include <stdio.h>
#include <string>
#include "Riostream.h"
#include "TVectorT.h"
#include "TString.h"
#include "TMath.h"

void MyEffWithDigis_all::Begin(TTree * /*tree*/)
{
  cout << " Begins " <<  endl;
 
}

void MyEffWithDigis_all::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(fmod(jentry,10000)==0)cout << " Processing entry : " << jentry <<" of "<<nentries<<endl;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     fChain->GetEntry(jentry); 
       
     // 4D Segments
     int NSegs=dtsegm4D_wheel->size();
     int NPhiSegs=0;
     int NZSegs=0;
     //if(fmod(jentry,10000)==0)cout << " Segments= " << NSegs << endl;
    
     // Loop in Segments to compute efficiencie L1Phi2 Different HV, Layer to Study
     for(int iseg=0;iseg<NSegs; iseg++)
       {
	 if((*dtsegm4D_hasPhi)[iseg] && (*dtsegm4D_phinhits)[iseg]>4 && fabs((*dtsegm4D_t0)[iseg])<50 )
	   {

	     int iwh=(*dtsegm4D_wheel)[iseg];
	     int is =(*dtsegm4D_sector)[iseg];
	     int ist=(*dtsegm4D_station)[iseg];
	     
	     TVectorF * VPhi_hitsSL = (TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(iseg); // Problems interactively
	     TVectorF * VPhi_hitsL = (TVectorF*)dtsegm4D_phi_hitsLayer->At(iseg); 
	     TVectorF * VPhi_hitsW = (TVectorF*)dtsegm4D_phi_hitsWire->At(iseg); 
	     TVectorF * VPhi_hitsPos = (TVectorF*)dtsegm4D_phi_hitsPos->At(iseg); 
	     float *  hitsSL_phi=VPhi_hitsSL->GetMatrixArray();
	     float *  hitsL_phi =VPhi_hitsL->GetMatrixArray();
	     float *  hitsW_phi =VPhi_hitsW->GetMatrixArray();
	     float *  hitsPos_phi =VPhi_hitsPos->GetMatrixArray();
	     int nhPhi  = (*dtsegm4D_phinhits)[iseg];
                   
	     // Fill Phihits in Segment
	     // hPhiHitsSeg[iwh+2][ist-1][is-1]->Fill(nhPhi); 
	
	     int totRelevantHits=0;
	     int WireFound=-1;int PosFound=-1;
	     int ReferencePhiHits=0; //All hits excluding the Phi2 L1 underStudy 
	     TVectorF * VLayersExpPos = (TVectorF*)dtsegm4D_hitsExpPos->At(iseg); // Problems interactively
	     TVectorF * VLayersExpWire = (TVectorF*)dtsegm4D_hitsExpWire->At(iseg); // Problems interactively
	     float * LayersExpPos=VLayersExpPos->GetMatrixArray();
	     float * LayersExpWire=VLayersExpWire->GetMatrixArray();
	     float residual=0;
	     

	     for(int ih=0;ih<nhPhi;ih++)
	       {// Loop in Hits
		 int isl=hitsSL_phi[ih];
		 int il= hitsL_phi[ih];

		 if(isl==3 && il==1){ // Phi2 Layer1
		   WireFound= hitsW_phi[ih];
		   PosFound= hitsPos_phi[ih];
		   residual=PosFound-LayersExpPos[8];
		 }
		 else ReferencePhiHits++;
	       }
	     if(ReferencePhiHits>4) //Hits in both SLs 
	       //if(ReferencePhiHits>6) //Hits in ALL layers of both SLs (except the layer under study) 
	       {
		 // All tracks
		 hAll[iwh+2][ist-1][is-1]->Fill(residual); 
	      
		 // Finding Hit 
		 if( WireFound ==  LayersExpWire[8] ) 
		   {
		     hFound[iwh+2][ist-1][is-1]->Fill(residual); 
		     hFoundDigi[iwh+2][ist-1][is-1]->Fill(residual); 
		   }
		 else
		   {// Look for possible unassociated hit inside the Digis
		     //cout<<"wire found != layer exp wire"<<endl;
		     int FoundDigi=0;
		     for(int nDigi=0; nDigi<Ndigis; nDigi++)
		       {
			 if( (*digi_wheel)[nDigi]== iwh         && (*digi_sector)[nDigi]==is    
			     && (*digi_station)[nDigi]==ist        && (*digi_sl)[nDigi]==3
			     && (*digi_layer)[nDigi]==1 && (*digi_wire)[nDigi] ==  LayersExpWire[8]
			     //&& (*digi_layer)[nDigi]==hitsL_phi[8] && abs((*digi_wire)[nDigi] - LayersExpWire[8])<2
			     )
			   FoundDigi++;
		       }
		     if(FoundDigi>0)  hFoundDigi[iwh+2][ist-1][is-1]->Fill(residual); 
		   }
	       }// End cut in number of hits Phi1+Phi2 excluding L1Phi2
	   }// End cut number of hits and T0
       }// End loop in Segments
   }// End loop in entries
}

void MyEffWithDigis_all::Terminate()
{

  //////////////////////////////////////////////////////////TABLE/////////////////////////////////////////////////////////////
  int nbins[5][4][14];
  float xntot[5][4][14];
  float xnfound[5][4][14];
  float xnfounddigi[5][4][14];
  ofstream efftxt (Form("eff_run%d.txt",runnumber));

  efftxt<<"    Wheel "<<"  "<<" MB "<<"  "<<" Sector "<<" "<< " Efficiency "<<""<<" Digi Efficiency "<<endl;
  efftxt<< " ................................................................................. " <<endl;

  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sector=0;sector<14;sector++){
	if(sector>11 && station!=3)continue;
	totalnumber[wheel][station][sector]=hAll[wheel][station][sector]  ->GetEntries();
	if(totalnumber[wheel][station][sector]==0){
 	  efftxt<<"missing chamber Wheel"<<wheel-2<<" MB"<<station+1<<" Sector"<<sector+1<<endl; 
 	  continue;
	}
	foundnumber[wheel][station][sector]=hFound[wheel][station][sector]->GetEntries();
	founddiginumber[wheel][station][sector]=hFoundDigi[wheel][station][sector]->GetEntries();
	Eff[wheel][station][sector]=100*foundnumber[wheel][station][sector]/totalnumber[wheel][station][sector];
	EffDigis[wheel][station][sector]=100*founddiginumber[wheel][station][sector]/totalnumber[wheel][station][sector];
	heff[station]->Fill(sector,wheel-3,EffDigis[wheel][station][sector]);

	efftxt<<"      "<<showpos<<wheel-2<<noshowpos<<"      "<<station+1<<"      "<<setfill('0') << setw(2)<< sector+1<<"      "<< fixed << setprecision(4)<<Eff[wheel][station][sector]<<"    "<<EffDigis[wheel][station][sector]<<endl;
      }
    }
  }
  cout << " Creating table   " <<Form("eff_run%d.txt",runnumber)  <<  endl;

 //////////////////////////////////////////////////////////COLOURED PLOTS/////////////////////////////////////////////////////////////
// 

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("5.3f");



  for(int i=0;i<3;i++) {
    heff[i]->GetYaxis()->SetNdivisions(5);
    heff[i]->GetXaxis()->SetNdivisions(12);  
    for(Int_t w=0;w<5;w++)heff[i]->GetYaxis()->SetBinLabel(w+1,wheel_label[w]);
    for(Int_t se=0;se<12;se++)heff[i]->GetXaxis()->SetBinLabel(se+1,sector_label[se]);
  }
  heff[3]->GetYaxis()->SetNdivisions(5);  
  heff[3]->GetXaxis()->SetNdivisions(14);
  for(Int_t w=0;w<5;w++)heff[3]->GetYaxis()->SetBinLabel(w+1,wheel_label[w]);
  for(Int_t se=0;se<14;se++)heff[3]->GetXaxis()->SetBinLabel(se+1,sector_label[se]);


  for(int i=0;i<4;i++) {
    heff[i]->SetTitle(Form("MB%d",i+1));
    heff[i]->GetXaxis()->SetTitle("Sector");
    heff[i]->GetYaxis()->SetTitle("Wheel");
    heff[i]->GetXaxis()->CenterTitle(1); 
    heff[i]->GetYaxis()->CenterTitle(1);
    heff[i]->GetXaxis()->SetTitleOffset(1.2);
    heff[i]->GetYaxis()->SetTitleOffset(1.4);
    heff[i]->SetStats(0);     
  }
  cout << " Creating canvases   "<<endl;

  TCanvas *can[4];
  for(int i=0;i<4;i++){

    can[i] = new TCanvas();
    can[i]->SetGrid();
    can[i]->SetLeftMargin(0.15);
    can[i]->SetBottomMargin(0.15);
    can[i]->SetRightMargin(0.25);
    can[i]->SetTopMargin(0.20);
    heff[i]->GetXaxis()->SetTitleFont(22);
    heff[i]->GetXaxis()->SetLabelFont(22);
    heff[i]->GetXaxis()->SetLabelSize(0.03);
    heff[i]->GetXaxis()->SetTitleSize(0.04);
    heff[i]->SetMaximum(100);
    heff[i]->SetMinimum(65);
    can[i]->cd();
    heff[i]->Draw("colz");
    heff[i]->Draw("TEXTsame45");
    can[i]->Print(Form("MB%d_efficiency.gif",i+1));
   }
  /////////////////////////////////////////////////EFFICIENCY VS PATH LENGTH/////////////////////////////////////
  
  // 4 Stations, 5 Wheels
  // Sector 4 only 
  double angle[5][4]={{38.6,44.2,49,52.8},{60,62.5,67.7,70.4},{90,90,90,90},{60,62.5,67.7,70.4},{38.6,44.2,49,52.8}};
  double path[5][4];
  
  // Calculate path lengths and fill with efficiencies
  for (int st=0; st<4;st++){
    int p=0;
    for (int w=0; w<5;w++){
      //if(w==2)continue;
      path[w][st]=TMath::Cos(angle[w][st]*TMath::Pi()/180)/TMath::Sin(angle[w][st]*TMath::Pi()/180);
      cout<<"Length of path for MB"<<st+1<<" Wheel"<<w-2<<" = "<<path[w][st]<<"  Efficiency = "<<EffDigis[w][st][3]<<endl;
      gr[st]->SetPoint(p,path[w][st],EffDigis[w][st][3]);
      gr[st]->SetMarkerColor(st+1);
      p++;
    }
  }
  
  TCanvas *c1= new TCanvas();
  for (int st=0; st<4;st++){mg->Add(gr[st]);gr[st]->SetMarkerSize(1.); gr[st]->SetMarkerStyle(20);}
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("<Path projection along wire> (effective height)");
  mg->GetYaxis()->SetTitle("Efficiency");
  mg->Fit("pol1");
  c1->SaveAs("effvspath.gif") ;
  c1->SaveAs("effvspath.root") ;
  
}

void MyEffWithDigis_all::SaveHistos(TFile * fhout)
{
  cout<<" Saving histograms"<<endl;
  
  for(int wheel=0;wheel<5;wheel++){
    for(int station=0;station<4;station++){
      for(int sector=0;sector<14;sector++){
	if(sector>11 && station!=3)continue;
	hAll[wheel][station][sector]->Write();
	hFound[wheel][station][sector]->Write();
	hFoundDigi[wheel][station][sector]->Write();
      }
    }
  }

for(int station=0;station<4;station++){
  heff[station]->Write();
 }



}
 
int main(int argc, char **argv){

 char infilename[400];
 char outfilename[400];

 if(argc < 3)
 {
    printf(" ***************************************************************************** \n" );
    printf(" ****    ERROR: Please enter the input directory and root filename \n" );
    printf(" ****    Usage: MyEffWithDigis_all.exe directorypath InputFilename.root \n");
    printf(" ****    EXITING PROGRAM!!!!! \n" );
    printf(" ***************************************************************************** \n" );
    exit(0);
 }
 else
 {
   sprintf(infilename,"%s/%s",argv[1],argv[2]);
   printf(" File to be used: %s \n",infilename);
 }

// TREE FROM INPUT FILE
   TChain chain("DTTree");
   chain.Add(infilename);

   MyEffWithDigis_all * EffMyAnalysis = new MyEffWithDigis_all();
   EffMyAnalysis->Begin(&chain);
   EffMyAnalysis->Init(&chain);

   Int_t nentries = (Int_t)chain.GetEntries();
   cout << " >>>>  " <<  nentries << " Entries to be processed"  <<  endl;
   EffMyAnalysis->Loop();

 // OUTPUT FILES
   TString inputfile =argv[2];
   TString outputfile=inputfile.ReplaceAll("DTTree_R","r");
   sprintf(outfilename,"eff_%s",outputfile.Data());
   TFile * fhistoout = new TFile (outfilename, "RECREATE");
   cout << " Creating root file " << outfilename  <<  endl;

   EffMyAnalysis->Terminate();
   EffMyAnalysis->SaveHistos(fhistoout);
   
   cout << " Ends " <<  endl;

} // END main







