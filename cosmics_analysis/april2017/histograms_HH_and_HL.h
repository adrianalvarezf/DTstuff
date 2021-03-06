

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  8 12:48:20 2016 by ROOT version 6.02/05
// from TTree DTTree/CMSSW DT tree
// found on file: /afs/cern.ch/user/a/alvareza/public/DTNtuple.root
//////////////////////////////////////////////////////////

#ifndef histograms_HH_and_HL_h
#define histograms_HH_and_HL_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "TClonesArray.h"
#include "vector"

#include <TH1F.h>
#include <TCanvas.h>

class histograms_HH_and_HL : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnumber;
   Int_t           lumiblock;
   Int_t           eventNumber;
   Float_t         timestamp;
   Int_t           bunchXing;
   Int_t           orbitNum;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_xxE;
   Float_t         PV_yyE;
   Float_t         PV_zzE;
   Float_t         PV_xyE;
   Float_t         PV_xzE;
   Float_t         PV_yzE;
   Float_t         PV_normchi2;
   Float_t         PV_Nvtx;
   Float_t         lumiperblock;
   Float_t         beam1Intensity;
   Float_t         beam2Intensity;

   int highhigh=0;
   int highlow=0;
   int worseq=0;
   int both=0;
   int totalsegm=0;

   TH1F *HH[5][4][14];  
   TH1F *HL[5][4][14]; 
   TH1F *div[5][4][14];
   TH1F *HH_fold[5][4][14];  
   TH1F *HL_fold[5][4][14];
   TH1F *div_fold[5][4][14]; 
  
   
   vector<TString> *hlt_path;
   vector<short>   *digi_wheel;
   vector<short>   *digi_sector;
   vector<short>   *digi_station;
   vector<short>   *digi_sl;
   vector<short>   *digi_layer;
   vector<short>   *digi_wire;
   vector<float>   *digi_time;
   vector<short>   *dtsegm4D_wheel;
   vector<short>   *dtsegm4D_sector;
   vector<short>   *dtsegm4D_station;
   vector<short>   *dtsegm4D_hasPhi;
   vector<short>   *dtsegm4D_hasZed;
   vector<float>   *dtsegm4D_x_pos_loc;
   vector<float>   *dtsegm4D_y_pos_loc;
   vector<float>   *dtsegm4D_z_pos_loc;
   vector<float>   *dtsegm4D_x_dir_loc;
   vector<float>   *dtsegm4D_y_dir_loc;
   vector<float>   *dtsegm4D_z_dir_loc;
   vector<float>   *dtsegm4D_cosx;
   vector<float>   *dtsegm4D_cosy;
   vector<float>   *dtsegm4D_cosz;
   vector<float>   *dtsegm4D_phi;
   vector<float>   *dtsegm4D_theta;
   vector<float>   *dtsegm4D_eta;
   vector<float>   *dtsegm4D_t0;
   vector<float>   *dtsegm4D_vdrift;
   vector<float>   *dtsegm4D_phinormchisq;
   vector<short>   *dtsegm4D_phinhits;
   vector<float>   *dtsegm4D_znormchisq;
   vector<short>   *dtsegm4D_znhits;

   TClonesArray    *dtsegm4D_hitsExpPos;
   TClonesArray    *dtsegm4D_hitsExpWire;
   TClonesArray    *dtsegm4D_phi_hitsPos;
   TClonesArray    *dtsegm4D_phi_hitsPosCh;
   TClonesArray    *dtsegm4D_phi_hitsPosErr;
   TClonesArray    *dtsegm4D_phi_hitsSide;
   TClonesArray    *dtsegm4D_phi_hitsWire;
   TClonesArray    *dtsegm4D_phi_hitsLayer;
   TClonesArray    *dtsegm4D_phi_hitsSuperLayer;
   TClonesArray    *dtsegm4D_z_hitsPos;
   TClonesArray    *dtsegm4D_z_hitsPosCh;
   TClonesArray    *dtsegm4D_z_hitsPosErr;
   TClonesArray    *dtsegm4D_z_hitsSide;
   TClonesArray    *dtsegm4D_z_hitsWire;
   TClonesArray    *dtsegm4D_z_hitsLayer;
   vector<short>   *cscsegm_ring;
   vector<short>   *cscsegm_chamber;
   vector<short>   *cscsegm_station;
   vector<float>   *cscsegm_cosx;
   vector<float>   *cscsegm_cosy;
   vector<float>   *cscsegm_cosz;
   vector<float>   *cscsegm_phi;
   vector<float>   *cscsegm_eta;
   vector<float>   *cscsegm_normchisq;
   vector<short>   *cscsegm_nRecHits;
   vector<short>   *ltTwinMuxIn_wheel;
   vector<short>   *ltTwinMuxIn_sector;
   vector<short>   *ltTwinMuxIn_station;
   vector<short>   *ltTwinMuxIn_quality;
   vector<short>   *ltTwinMuxIn_bx;
   vector<short>   *ltTwinMuxIn_phi;
   vector<short>   *ltTwinMuxIn_phiB;
   vector<short>   *ltTwinMuxIn_is2nd;
   vector<short>   *ltTwinMux_thWheel;
   vector<short>   *ltTwinMux_thSector;
   vector<short>   *ltTwinMux_thStation;
   vector<short>   *ltTwinMux_thBx;
   vector<short>   *ltTwinMux_thHits;
   vector<short>   *Mu_isMuGlobal;
   vector<short>   *Mu_isMuTracker;
   vector<int>     *Mu_numberOfChambers_sta;
   vector<int>     *Mu_numberOfMatches_sta;
   vector<int>     *Mu_numberOfHits_sta;
   vector<int>     *Mu_segmentIndex_sta;
   vector<float>   *Mu_px;
   vector<float>   *Mu_py;
   vector<float>   *Mu_pz;
   vector<float>   *Mu_phi;
   vector<float>   *Mu_eta;
   vector<short>   *Mu_recHitsSize;
   vector<float>   *Mu_normchi2_sta;
   vector<short>   *Mu_charge;
   vector<float>   *Mu_dxy_sta;
   vector<float>   *Mu_dz_sta;
   vector<float>   *Mu_normchi2_glb;
   vector<float>   *Mu_dxy_glb;
   vector<float>   *Mu_dz_glb;
   vector<int>     *Mu_numberOfPixelHits_glb;
   vector<int>     *Mu_numberOfTrackerHits_glb;
   vector<float>   *Mu_tkIsoR03_glb;
   vector<float>   *Mu_ntkIsoR03_glb;
   vector<float>   *Mu_emIsoR03_glb;
   vector<float>   *Mu_hadIsoR03_glb;
   vector<float>   *STAMu_caloCompatibility;
   vector<float>   *Mu_z_mb2_mu;
   vector<float>   *Mu_phi_mb2_mu;
   vector<float>   *Mu_pseta_mb2_mu;
   vector<short>   *gmt_bx;
   vector<float>   *gmt_phi;
   vector<float>   *gmt_eta;
   vector<float>   *gmt_pt;
   vector<short>   *gmt_qual;
   vector<short>   *gmt_detector;
   vector<short>   *gmt_cands_fwd;
   vector<short>   *gmt_cands_isRpc;
   vector<short>   *gmt_cands_bx;
   vector<float>   *gmt_cands_phi;
   vector<float>   *gmt_cands_eta;
   vector<float>   *gmt_cands_pt;
   vector<short>   *gmt_cands_qual;
   vector<short>   *gmt_cands_ismatched;
   vector<short>   *gt_algo_bit;
   vector<short>   *gt_algo_bx;
   vector<short>   *gt_tt_bit;
   vector<short>   *gt_tt_bx;
   vector<int>     *rpc_region;
   vector<int>     *rpc_clusterSize;
   vector<int>     *rpc_strip;
   vector<int>     *rpc_bx;
   vector<int>     *rpc_station;
   vector<int>     *rpc_sector;
   vector<int>     *rpc_layer;
   vector<int>     *rpc_subsector;
   vector<int>     *rpc_roll;
   vector<int>     *rpc_ring;
   Short_t         Ndigis;
   Short_t         Ndtsegments;
   Short_t         Ncscsegments;
   Short_t         Nmuons;
   Short_t         Ngmt;
   Short_t         Ngmtcands;
   Short_t         Ngtalgo;
   Short_t         Ngttechtrig;
   Short_t         Nhlt;
   Short_t         NrpcRecHits;

   // List of branches
   TBranch        *b_runnumber;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_bunchXing;   //!
   TBranch        *b_orbitNum;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_xxE;   //!
   TBranch        *b_PV_yyE;   //!
   TBranch        *b_PV_zzE;   //!
   TBranch        *b_PV_xyE;   //!
   TBranch        *b_PV_xzE;   //!
   TBranch        *b_PV_yzE;   //!
   TBranch        *b_PV_normch2;   //!
   TBranch        *b_PV_Nvtx;   //!
   TBranch        *b_lumiperblock;   //!
   TBranch        *b_beam1Intensity;   //!
   TBranch        *b_beam2Intensity;   //!
   TBranch        *b_hlt_path;   //!
   TBranch        *b_digi_wheel;   //!
   TBranch        *b_digi_sector;   //!
   TBranch        *b_digi_station;   //!
   TBranch        *b_digi_sl;   //!
   TBranch        *b_digi_layer;   //!
   TBranch        *b_digi_wire;   //!
   TBranch        *b_digi_time;   //!
   TBranch        *b_dtsegm4D_wheel;   //!
   TBranch        *b_dtsegm4D_sector;   //!
   TBranch        *b_dtsegm4D_station;   //!
   TBranch        *b_dtsegm4D_hasPhi;   //!
   TBranch        *b_dtsegm4D_hasZed;   //!
   TBranch        *b_dtsegm4D_x_pos_loc;   //!
   TBranch        *b_dtsegm4D_y_pos_loc;   //!
   TBranch        *b_dtsegm4D_z_pos_loc;   //!
   TBranch        *b_dtsegm4D_x_dir_loc;   //!
   TBranch        *b_dtsegm4D_y_dir_loc;   //!
   TBranch        *b_dtsegm4D_z_dir_loc;   //!
   TBranch        *b_dtsegm4D_cosx;   //!
   TBranch        *b_dtsegm4D_cosy;   //!
   TBranch        *b_dtsegm4D_cosz;   //!
   TBranch        *b_dtsegm4D_phi;   //!
   TBranch        *b_dtsegm4D_theta;   //!
   TBranch        *b_dtsegm4D_eta;   //!
   TBranch        *b_dtsegm4D_t0;   //!
   TBranch        *b_dtsegm4D_vdrift;   //!
   TBranch        *b_dtsegm4D_phinormchisq;   //!
   TBranch        *b_dtsegm4D_phinhits;   //!
   TBranch        *b_dtsegm4D_znormchisq;   //!
   TBranch        *b_dtsegm4D_znhits;   //!
   TBranch        *b_dtsegm4D_hitsExpPos;   //!
   TBranch        *b_dtsegm4D_hitsExpWire;   //!
   TBranch        *b_dtsegm4D_phi_hitsPos;   //!
   TBranch        *b_dtsegm4D_phi_hitsPosCh;   //!
   TBranch        *b_dtsegm4D_phi_hitsPosErr;   //!
   TBranch        *b_dtsegm4D_phi_hitsSide;   //!
   TBranch        *b_dtsegm4D_phi_hitsWire;   //!
   TBranch        *b_dtsegm4D_phi_hitsLayer;   //!
   TBranch        *b_dtsegm4D_phi_hitsSuperLayer;   //!
   TBranch        *b_dtsegm4D_z_hitsPos;   //!
   TBranch        *b_dtsegm4D_z_hitsPosCh;   //!
   TBranch        *b_dtsegm4D_z_hitsPosErr;   //!
   TBranch        *b_dtsegm4D_z_hitsSide;   //!
   TBranch        *b_dtsegm4D_z_hitsWire;   //!
   TBranch        *b_dtsegm4D_z_hitsLayer;   //!
   TBranch        *b_cscsegm_ring;   //!
   TBranch        *b_cscsegm_chamber;   //!
   TBranch        *b_cscsegm_station;   //!
   TBranch        *b_cscsegm_cosx;   //!
   TBranch        *b_cscsegm_cosy;   //!
   TBranch        *b_cscsegm_cosz;   //!
   TBranch        *b_cscsegm_phi;   //!
   TBranch        *b_cscsegm_eta;   //!
   TBranch        *b_cscsegm_normchisq;   //!
   TBranch        *b_cscsegm_nRecHits;   //!
   TBranch        *b_ltTwinMuxIn_wheel;//!
   TBranch        *b_ltTwinMuxIn_sector;//!
   TBranch        *b_ltTwinMuxIn_station;//!
   TBranch        *b_ltTwinMuxIn_quality;//!
   TBranch        *b_ltTwinMuxIn_bx;//!
   TBranch        *b_ltTwinMuxIn_phi;//!
   TBranch        *b_ltTwinMuxIn_phiB;//!
   TBranch        *b_ltTwinMuxIn_is2nd;//!
   TBranch        *b_ltTwinMux_thWheel;//!
   TBranch        *b_ltTwinMux_thSector;//!
   TBranch        *b_ltTwinMux_thStation;//!
   TBranch        *b_ltTwinMux_thBx;//!
   TBranch        *b_ltTwinMux_thHits;//!
   TBranch        *b_Mu_isMuGlobal;   //!
   TBranch        *b_Mu_isMuTracker;   //!
   TBranch        *b_Mu_numberOfChambers_sta;   //!
   TBranch        *b_Mu_numberOfMatches_sta;   //!
   TBranch        *b_Mu_numberOfHits_sta;   //!
   TBranch        *b_Mu_segmentIndex_sta;   //!
   TBranch        *b_Mu_px;   //!
   TBranch        *b_Mu_py;   //!
   TBranch        *b_Mu_pz;   //!
   TBranch        *b_Mu_phi;   //!
   TBranch        *b_Mu_eta;   //!
   TBranch        *b_Mu_recHitsSize;   //!
   TBranch        *b_Mu_normchi2_sta;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_Mu_dxy_sta;   //!
   TBranch        *b_Mu_dz_sta;   //!
   TBranch        *b_Mu_normchi2_glb;   //!
   TBranch        *b_Mu_dxy_glb;   //!
   TBranch        *b_Mu_dz_glb;   //!
   TBranch        *b_Mu_numberOfPixelHits_glb;   //!
   TBranch        *b_Mu_numberOfTrackerHits_glb;   //!
   TBranch        *b_Mu_tkIsoR03_glb;   //!
   TBranch        *b_Mu_ntkIsoR03_glb;   //!
   TBranch        *b_Mu_emIsoR03_glb;   //!
   TBranch        *b_Mu_hadIsoR03_glb;   //!
   TBranch        *b_STAMu_caloCompatibility;   //!
   TBranch        *b_Mu_z_mb2_mu;   //!
   TBranch        *b_Mu_phi_mb2_mu;   //!
   TBranch        *b_Mu_pseta_mb2_mu;   //!
   TBranch        *b_gmt_bx;   //!
   TBranch        *b_gmt_phi;   //!
   TBranch        *b_gmt_eta;   //!
   TBranch        *b_gmt_pt;   //!
   TBranch        *b_gmt_qual;   //!
   TBranch        *b_gmt_detector;   //!
   TBranch        *b_gmt_cands_fwd;   //!
   TBranch        *b_gmt_cands_isRpc;   //!
   TBranch        *b_gmt_cands_bx;   //!
   TBranch        *b_gmt_cands_phi;   //!
   TBranch        *b_gmt_cands_eta;   //!
   TBranch        *b_gmt_cands_pt;   //!
   TBranch        *b_gmt_cands_qual;   //!
   TBranch        *b_gmt_cands_ismatched;   //!
   TBranch        *b_gt_algo_bit;   //!
   TBranch        *b_gt_algo_bx;   //!
   TBranch        *b_gt_tt_bit;   //!
   TBranch        *b_gt_tt_bx;   //!
   TBranch        *b_rpc_region;   //!
   TBranch        *b_rpc_clusterSize;   //!
   TBranch        *b_rpc_strip;   //!
   TBranch        *b_rpc_bx;   //!
   TBranch        *b_rpc_station;   //!
   TBranch        *b_rpc_sector;   //!
   TBranch        *b_rpc_layer;   //!
   TBranch        *b_rpc_subsector;   //!
   TBranch        *b_rpc_roll;   //!
   TBranch        *b_rpc_ring;   //!
   TBranch        *b_Ndigis;   //!
   TBranch        *b_Ndtsegments;   //!
   TBranch        *b_Ncscsegments;   //!
   TBranch        *b_Nmuons;   //!
   TBranch        *b_Ngmt;   //!
   TBranch        *b_Ngmtcands;   //!
   TBranch        *b_Ngtalgo;   //!
   TBranch        *b_Ngttt;   //!
   TBranch        *b_Nhlt;   //!
   TBranch        *b_NrpcRecHits;   //!

   histograms_HH_and_HL(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~histograms_HH_and_HL() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    Terminate();

   ClassDef(histograms_HH_and_HL,0);
};

#endif

#ifdef histograms_HH_and_HL_cxx
void histograms_HH_and_HL::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hlt_path = 0;
   digi_wheel = 0;
   digi_sector = 0;
   digi_station = 0;
   digi_sl = 0;
   digi_layer = 0;
   digi_wire = 0;
   digi_time = 0;
   dtsegm4D_wheel = 0;
   dtsegm4D_sector = 0;
   dtsegm4D_station = 0;
   dtsegm4D_hasPhi = 0;
   dtsegm4D_hasZed = 0;
   dtsegm4D_x_pos_loc = 0;
   dtsegm4D_y_pos_loc = 0;
   dtsegm4D_z_pos_loc = 0;
   dtsegm4D_x_dir_loc = 0;
   dtsegm4D_y_dir_loc = 0;
   dtsegm4D_z_dir_loc = 0;
   dtsegm4D_cosx = 0;
   dtsegm4D_cosy = 0;
   dtsegm4D_cosz = 0;
   dtsegm4D_phi = 0;
   dtsegm4D_theta = 0;
   dtsegm4D_eta = 0;
   dtsegm4D_t0 = 0;
   dtsegm4D_vdrift = 0;
   dtsegm4D_phinormchisq = 0;
   dtsegm4D_phinhits = 0;
   dtsegm4D_znormchisq = 0;
   dtsegm4D_znhits = 0;
   dtsegm4D_hitsExpPos = 0;
   dtsegm4D_hitsExpWire = 0;
   dtsegm4D_phi_hitsPos = 0;
   dtsegm4D_phi_hitsPosCh = 0;
   dtsegm4D_phi_hitsPosErr = 0;
   dtsegm4D_phi_hitsSide = 0;
   dtsegm4D_phi_hitsWire = 0;
   dtsegm4D_phi_hitsLayer = 0;
   dtsegm4D_phi_hitsSuperLayer = 0;
   dtsegm4D_z_hitsPos = 0;
   dtsegm4D_z_hitsPosCh = 0;
   dtsegm4D_z_hitsPosErr = 0;
   dtsegm4D_z_hitsSide = 0;
   dtsegm4D_z_hitsWire = 0;
   dtsegm4D_z_hitsLayer = 0;
   cscsegm_ring = 0;
   cscsegm_chamber = 0;
   cscsegm_station = 0;
   cscsegm_cosx = 0;
   cscsegm_cosy = 0;
   cscsegm_cosz = 0;
   cscsegm_phi = 0;
   cscsegm_eta = 0;
   cscsegm_normchisq = 0;
   cscsegm_nRecHits = 0;
   ltTwinMuxIn_wheel = 0;
   ltTwinMuxIn_sector = 0;
   ltTwinMuxIn_station = 0;
   ltTwinMuxIn_quality = 0;
   ltTwinMuxIn_bx = 0;
   ltTwinMuxIn_phi = 0;
   ltTwinMuxIn_phiB = 0;
   ltTwinMuxIn_is2nd = 0;
   ltTwinMux_thWheel = 0;
   ltTwinMux_thSector = 0;
   ltTwinMux_thStation = 0;
   ltTwinMux_thBx = 0;
   ltTwinMux_thHits = 0;
   Mu_isMuGlobal = 0;
   Mu_isMuTracker = 0;
   Mu_numberOfChambers_sta = 0;
   Mu_numberOfMatches_sta = 0;
   Mu_numberOfHits_sta = 0;
   Mu_segmentIndex_sta = 0;
   Mu_px = 0;
   Mu_py = 0;
   Mu_pz = 0;
   Mu_phi = 0;
   Mu_eta = 0;
   Mu_recHitsSize = 0;
   Mu_normchi2_sta = 0;
   Mu_charge = 0;
   Mu_dxy_sta = 0;
   Mu_dz_sta = 0;
   Mu_normchi2_glb = 0;
   Mu_dxy_glb = 0;
   Mu_dz_glb = 0;
   Mu_numberOfPixelHits_glb = 0;
   Mu_numberOfTrackerHits_glb = 0;
   Mu_tkIsoR03_glb = 0;
   Mu_ntkIsoR03_glb = 0;
   Mu_emIsoR03_glb = 0;
   Mu_hadIsoR03_glb = 0;
   STAMu_caloCompatibility = 0;
   Mu_z_mb2_mu = 0;
   Mu_phi_mb2_mu = 0;
   Mu_pseta_mb2_mu = 0;
   gmt_bx = 0;
   gmt_phi = 0;
   gmt_eta = 0;
   gmt_pt = 0;
   gmt_qual = 0;
   gmt_detector = 0;
   gmt_cands_fwd = 0;
   gmt_cands_isRpc = 0;
   gmt_cands_bx = 0;
   gmt_cands_phi = 0;
   gmt_cands_eta = 0;
   gmt_cands_pt = 0;
   gmt_cands_qual = 0;
   gmt_cands_ismatched = 0;
   gt_algo_bit = 0;
   gt_algo_bx = 0;
   gt_tt_bit = 0;
   gt_tt_bx = 0;
   rpc_region = 0;
   rpc_clusterSize = 0;
   rpc_strip = 0;
   rpc_bx = 0;
   rpc_station = 0;
   rpc_sector = 0;
   rpc_layer = 0;
   rpc_subsector = 0;
   rpc_roll = 0;
   rpc_ring = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("bunchXing", &bunchXing, &b_bunchXing);
   fChain->SetBranchAddress("orbitNum", &orbitNum, &b_orbitNum);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_xxE", &PV_xxE, &b_PV_xxE);
   fChain->SetBranchAddress("PV_yyE", &PV_yyE, &b_PV_yyE);
   fChain->SetBranchAddress("PV_zzE", &PV_zzE, &b_PV_zzE);
   fChain->SetBranchAddress("PV_xyE", &PV_xyE, &b_PV_xyE);
   fChain->SetBranchAddress("PV_xzE", &PV_xzE, &b_PV_xzE);
   fChain->SetBranchAddress("PV_yzE", &PV_yzE, &b_PV_yzE);
   fChain->SetBranchAddress("PV_normchi2", &PV_normchi2, &b_PV_normch2);
   fChain->SetBranchAddress("PV_Nvtx", &PV_Nvtx, &b_PV_Nvtx);
   fChain->SetBranchAddress("lumiperblock", &lumiperblock, &b_lumiperblock);
   fChain->SetBranchAddress("beam1Intensity", &beam1Intensity, &b_beam1Intensity);
   fChain->SetBranchAddress("beam2Intensity", &beam2Intensity, &b_beam2Intensity);
   fChain->SetBranchAddress("hlt_path", &hlt_path, &b_hlt_path);
   fChain->SetBranchAddress("digi_wheel", &digi_wheel, &b_digi_wheel);
   fChain->SetBranchAddress("digi_sector", &digi_sector, &b_digi_sector);
   fChain->SetBranchAddress("digi_station", &digi_station, &b_digi_station);
   fChain->SetBranchAddress("digi_sl", &digi_sl, &b_digi_sl);
   fChain->SetBranchAddress("digi_layer", &digi_layer, &b_digi_layer);
   fChain->SetBranchAddress("digi_wire", &digi_wire, &b_digi_wire);
   fChain->SetBranchAddress("digi_time", &digi_time, &b_digi_time);
   fChain->SetBranchAddress("dtsegm4D_wheel", &dtsegm4D_wheel, &b_dtsegm4D_wheel);
   fChain->SetBranchAddress("dtsegm4D_sector", &dtsegm4D_sector, &b_dtsegm4D_sector);
   fChain->SetBranchAddress("dtsegm4D_station", &dtsegm4D_station, &b_dtsegm4D_station);
   fChain->SetBranchAddress("dtsegm4D_hasPhi", &dtsegm4D_hasPhi, &b_dtsegm4D_hasPhi);
   fChain->SetBranchAddress("dtsegm4D_hasZed", &dtsegm4D_hasZed, &b_dtsegm4D_hasZed);
   fChain->SetBranchAddress("dtsegm4D_x_pos_loc", &dtsegm4D_x_pos_loc, &b_dtsegm4D_x_pos_loc);
   fChain->SetBranchAddress("dtsegm4D_y_pos_loc", &dtsegm4D_y_pos_loc, &b_dtsegm4D_y_pos_loc);
   fChain->SetBranchAddress("dtsegm4D_z_pos_loc", &dtsegm4D_z_pos_loc, &b_dtsegm4D_z_pos_loc);
   fChain->SetBranchAddress("dtsegm4D_x_dir_loc", &dtsegm4D_x_dir_loc, &b_dtsegm4D_x_dir_loc);
   fChain->SetBranchAddress("dtsegm4D_y_dir_loc", &dtsegm4D_y_dir_loc, &b_dtsegm4D_y_dir_loc);
   fChain->SetBranchAddress("dtsegm4D_z_dir_loc", &dtsegm4D_z_dir_loc, &b_dtsegm4D_z_dir_loc);
   fChain->SetBranchAddress("dtsegm4D_cosx", &dtsegm4D_cosx, &b_dtsegm4D_cosx);
   fChain->SetBranchAddress("dtsegm4D_cosy", &dtsegm4D_cosy, &b_dtsegm4D_cosy);
   fChain->SetBranchAddress("dtsegm4D_cosz", &dtsegm4D_cosz, &b_dtsegm4D_cosz);
   fChain->SetBranchAddress("dtsegm4D_phi", &dtsegm4D_phi, &b_dtsegm4D_phi);
   fChain->SetBranchAddress("dtsegm4D_theta", &dtsegm4D_theta, &b_dtsegm4D_theta);
   fChain->SetBranchAddress("dtsegm4D_eta", &dtsegm4D_eta, &b_dtsegm4D_eta);
   fChain->SetBranchAddress("dtsegm4D_t0", &dtsegm4D_t0, &b_dtsegm4D_t0);
   fChain->SetBranchAddress("dtsegm4D_vdrift", &dtsegm4D_vdrift, &b_dtsegm4D_vdrift);
   fChain->SetBranchAddress("dtsegm4D_phinormchisq", &dtsegm4D_phinormchisq, &b_dtsegm4D_phinormchisq);
   fChain->SetBranchAddress("dtsegm4D_phinhits", &dtsegm4D_phinhits, &b_dtsegm4D_phinhits);
   fChain->SetBranchAddress("dtsegm4D_znormchisq", &dtsegm4D_znormchisq, &b_dtsegm4D_znormchisq);
   fChain->SetBranchAddress("dtsegm4D_znhits", &dtsegm4D_znhits, &b_dtsegm4D_znhits);
   fChain->SetBranchAddress("dtsegm4D_hitsExpPos", &dtsegm4D_hitsExpPos, &b_dtsegm4D_hitsExpPos);
   fChain->SetBranchAddress("dtsegm4D_hitsExpWire",  &dtsegm4D_hitsExpWire, &b_dtsegm4D_hitsExpWire);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsPos", &dtsegm4D_phi_hitsPos, &b_dtsegm4D_phi_hitsPos);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsPosCh", &dtsegm4D_phi_hitsPosCh, &b_dtsegm4D_phi_hitsPosCh);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsPosErr", &dtsegm4D_phi_hitsPosErr, &b_dtsegm4D_phi_hitsPosErr);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsSide", &dtsegm4D_phi_hitsSide, &b_dtsegm4D_phi_hitsSide);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsWire", &dtsegm4D_phi_hitsWire, &b_dtsegm4D_phi_hitsWire);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsLayer", &dtsegm4D_phi_hitsLayer, &b_dtsegm4D_phi_hitsLayer);
   fChain->SetBranchAddress("dtsegm4D_phi_hitsSuperLayer", &dtsegm4D_phi_hitsSuperLayer, &b_dtsegm4D_phi_hitsSuperLayer);
   fChain->SetBranchAddress("dtsegm4D_z_hitsPos", &dtsegm4D_z_hitsPos, &b_dtsegm4D_z_hitsPos);
   fChain->SetBranchAddress("dtsegm4D_z_hitsPosCh", &dtsegm4D_z_hitsPosCh, &b_dtsegm4D_z_hitsPosCh);
   fChain->SetBranchAddress("dtsegm4D_z_hitsPosErr", &dtsegm4D_z_hitsPosErr, &b_dtsegm4D_z_hitsPosErr);
   fChain->SetBranchAddress("dtsegm4D_z_hitsSide", &dtsegm4D_z_hitsSide, &b_dtsegm4D_z_hitsSide);
   fChain->SetBranchAddress("dtsegm4D_z_hitsWire", &dtsegm4D_z_hitsWire, &b_dtsegm4D_z_hitsWire);
   fChain->SetBranchAddress("dtsegm4D_z_hitsLayer", &dtsegm4D_z_hitsLayer, &b_dtsegm4D_z_hitsLayer);
   fChain->SetBranchAddress("cscsegm_ring", &cscsegm_ring, &b_cscsegm_ring);
   fChain->SetBranchAddress("cscsegm_chamber", &cscsegm_chamber, &b_cscsegm_chamber);
   fChain->SetBranchAddress("cscsegm_station", &cscsegm_station, &b_cscsegm_station);
   fChain->SetBranchAddress("cscsegm_cosx", &cscsegm_cosx, &b_cscsegm_cosx);
   fChain->SetBranchAddress("cscsegm_cosy", &cscsegm_cosy, &b_cscsegm_cosy);
   fChain->SetBranchAddress("cscsegm_cosz", &cscsegm_cosz, &b_cscsegm_cosz);
   fChain->SetBranchAddress("cscsegm_phi", &cscsegm_phi, &b_cscsegm_phi);
   fChain->SetBranchAddress("cscsegm_eta", &cscsegm_eta, &b_cscsegm_eta);
   fChain->SetBranchAddress("cscsegm_normchisq", &cscsegm_normchisq, &b_cscsegm_normchisq);
   fChain->SetBranchAddress("cscsegm_nRecHits", &cscsegm_nRecHits, &b_cscsegm_nRecHits);
   fChain->SetBranchAddress("ltTwinMuxIn_wheel",&ltTwinMuxIn_wheel,&b_ltTwinMuxIn_wheel);
   fChain->SetBranchAddress("ltTwinMuxIn_sector",&ltTwinMuxIn_sector,&b_ltTwinMuxIn_sector);
   fChain->SetBranchAddress("ltTwinMuxIn_station",&ltTwinMuxIn_station,&b_ltTwinMuxIn_station);
   fChain->SetBranchAddress("ltTwinMuxIn_quality",&ltTwinMuxIn_quality,&b_ltTwinMuxIn_quality);
   fChain->SetBranchAddress("ltTwinMuxIn_bx",&ltTwinMuxIn_bx,&b_ltTwinMuxIn_bx);
   fChain->SetBranchAddress("ltTwinMuxIn_phi",&ltTwinMuxIn_phi,&b_ltTwinMuxIn_phi);
   fChain->SetBranchAddress("ltTwinMuxIn_phiB",&ltTwinMuxIn_phiB,&b_ltTwinMuxIn_phiB);
   fChain->SetBranchAddress("ltTwinMuxIn_is2nd",&ltTwinMuxIn_is2nd,&b_ltTwinMuxIn_is2nd);
   fChain->SetBranchAddress("ltTwinMux_thWheel",&ltTwinMux_thWheel,&b_ltTwinMux_thWheel);
   fChain->SetBranchAddress("ltTwinMux_thSector",&ltTwinMux_thSector,&b_ltTwinMux_thSector);
   fChain->SetBranchAddress("ltTwinMux_thStation",&ltTwinMux_thStation,&b_ltTwinMux_thStation);
   fChain->SetBranchAddress("ltTwinMux_thBx",&ltTwinMux_thBx,&b_ltTwinMux_thBx);
   fChain->SetBranchAddress("ltTwinMux_thHits",&ltTwinMux_thHits,&b_ltTwinMux_thHits);
   fChain->SetBranchAddress("Mu_isMuGlobal", &Mu_isMuGlobal, &b_Mu_isMuGlobal);
   fChain->SetBranchAddress("Mu_isMuTracker", &Mu_isMuTracker, &b_Mu_isMuTracker);
   fChain->SetBranchAddress("Mu_numberOfChambers_sta", &Mu_numberOfChambers_sta, &b_Mu_numberOfChambers_sta);
   fChain->SetBranchAddress("Mu_numberOfMatches_sta", &Mu_numberOfMatches_sta, &b_Mu_numberOfMatches_sta);
   fChain->SetBranchAddress("Mu_numberOfHits_sta", &Mu_numberOfHits_sta, &b_Mu_numberOfHits_sta);
   fChain->SetBranchAddress("Mu_segmentIndex_sta", &Mu_segmentIndex_sta, &b_Mu_segmentIndex_sta);
   fChain->SetBranchAddress("Mu_px", &Mu_px, &b_Mu_px);
   fChain->SetBranchAddress("Mu_py", &Mu_py, &b_Mu_py);
   fChain->SetBranchAddress("Mu_pz", &Mu_pz, &b_Mu_pz);
   fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
   fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
   fChain->SetBranchAddress("Mu_recHitsSize", &Mu_recHitsSize, &b_Mu_recHitsSize);
   fChain->SetBranchAddress("Mu_normchi2_sta", &Mu_normchi2_sta, &b_Mu_normchi2_sta);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("Mu_dxy_sta", &Mu_dxy_sta, &b_Mu_dxy_sta);
   fChain->SetBranchAddress("Mu_dz_sta", &Mu_dz_sta, &b_Mu_dz_sta);
   fChain->SetBranchAddress("Mu_normchi2_glb", &Mu_normchi2_glb, &b_Mu_normchi2_glb);
   fChain->SetBranchAddress("Mu_dxy_glb", &Mu_dxy_glb, &b_Mu_dxy_glb);
   fChain->SetBranchAddress("Mu_dz_glb", &Mu_dz_glb, &b_Mu_dz_glb);
   fChain->SetBranchAddress("Mu_numberOfPixelHits_glb", &Mu_numberOfPixelHits_glb, &b_Mu_numberOfPixelHits_glb);
   fChain->SetBranchAddress("Mu_numberOfTrackerHits_glb", &Mu_numberOfTrackerHits_glb, &b_Mu_numberOfTrackerHits_glb);
   fChain->SetBranchAddress("Mu_tkIsoR03_glb", &Mu_tkIsoR03_glb, &b_Mu_tkIsoR03_glb);
   fChain->SetBranchAddress("Mu_ntkIsoR03_glb", &Mu_ntkIsoR03_glb, &b_Mu_ntkIsoR03_glb);
   fChain->SetBranchAddress("Mu_emIsoR03_glb", &Mu_emIsoR03_glb, &b_Mu_emIsoR03_glb);
   fChain->SetBranchAddress("Mu_hadIsoR03_glb", &Mu_hadIsoR03_glb, &b_Mu_hadIsoR03_glb);
   fChain->SetBranchAddress("STAMu_caloCompatibility", &STAMu_caloCompatibility, &b_STAMu_caloCompatibility);
   fChain->SetBranchAddress("Mu_z_mb2_mu", &Mu_z_mb2_mu, &b_Mu_z_mb2_mu);
   fChain->SetBranchAddress("Mu_phi_mb2_mu", &Mu_phi_mb2_mu, &b_Mu_phi_mb2_mu);
   fChain->SetBranchAddress("Mu_pseta_mb2_mu", &Mu_pseta_mb2_mu, &b_Mu_pseta_mb2_mu);
   fChain->SetBranchAddress("gmt_bx", &gmt_bx, &b_gmt_bx);
   fChain->SetBranchAddress("gmt_phi", &gmt_phi, &b_gmt_phi);
   fChain->SetBranchAddress("gmt_eta", &gmt_eta, &b_gmt_eta);
   fChain->SetBranchAddress("gmt_pt", &gmt_pt, &b_gmt_pt);
   fChain->SetBranchAddress("gmt_qual", &gmt_qual, &b_gmt_qual);
   fChain->SetBranchAddress("gmt_detector", &gmt_detector, &b_gmt_detector);
   fChain->SetBranchAddress("gmt_cands_fwd", &gmt_cands_fwd, &b_gmt_cands_fwd);
   fChain->SetBranchAddress("gmt_cands_isRpc", &gmt_cands_isRpc, &b_gmt_cands_isRpc);
   fChain->SetBranchAddress("gmt_cands_bx", &gmt_cands_bx, &b_gmt_cands_bx);
   fChain->SetBranchAddress("gmt_cands_phi", &gmt_cands_phi, &b_gmt_cands_phi);
   fChain->SetBranchAddress("gmt_cands_eta", &gmt_cands_eta, &b_gmt_cands_eta);
   fChain->SetBranchAddress("gmt_cands_pt", &gmt_cands_pt, &b_gmt_cands_pt);
   fChain->SetBranchAddress("gmt_cands_qual", &gmt_cands_qual, &b_gmt_cands_qual);
   fChain->SetBranchAddress("gmt_cands_ismatched", &gmt_cands_ismatched, &b_gmt_cands_ismatched);
   fChain->SetBranchAddress("gt_algo_bit", &gt_algo_bit, &b_gt_algo_bit);
   fChain->SetBranchAddress("gt_algo_bx", &gt_algo_bx, &b_gt_algo_bx);
   fChain->SetBranchAddress("gt_tt_bit", &gt_tt_bit, &b_gt_tt_bit);
   fChain->SetBranchAddress("gt_tt_bx", &gt_tt_bx, &b_gt_tt_bx);
   fChain->SetBranchAddress("rpc_region", &rpc_region, &b_rpc_region);
   fChain->SetBranchAddress("rpc_clusterSize", &rpc_clusterSize, &b_rpc_clusterSize);
   fChain->SetBranchAddress("rpc_strip", &rpc_strip, &b_rpc_strip);
   fChain->SetBranchAddress("rpc_bx", &rpc_bx, &b_rpc_bx);
   fChain->SetBranchAddress("rpc_station", &rpc_station, &b_rpc_station);
   fChain->SetBranchAddress("rpc_sector", &rpc_sector, &b_rpc_sector);
   fChain->SetBranchAddress("rpc_layer", &rpc_layer, &b_rpc_layer);
   fChain->SetBranchAddress("rpc_subsector", &rpc_subsector, &b_rpc_subsector);
   fChain->SetBranchAddress("rpc_roll", &rpc_roll, &b_rpc_roll);
   fChain->SetBranchAddress("rpc_ring", &rpc_ring, &b_rpc_ring);
   fChain->SetBranchAddress("Ndigis", &Ndigis, &b_Ndigis);
   fChain->SetBranchAddress("Ndtsegments", &Ndtsegments, &b_Ndtsegments);
   fChain->SetBranchAddress("Ncscsegments", &Ncscsegments, &b_Ncscsegments);
   fChain->SetBranchAddress("Nmuons", &Nmuons, &b_Nmuons);
   fChain->SetBranchAddress("Ngmt", &Ngmt, &b_Ngmt);
   fChain->SetBranchAddress("Ngmtcands", &Ngmtcands, &b_Ngmtcands);
   fChain->SetBranchAddress("Ngtalgo", &Ngtalgo, &b_Ngtalgo);
   fChain->SetBranchAddress("Ngttechtrig", &Ngttechtrig, &b_Ngttt);
   fChain->SetBranchAddress("Nhlt", &Nhlt, &b_Nhlt);
   fChain->SetBranchAddress("NrpcRecHits", &NrpcRecHits, &b_NrpcRecHits);
   
   for(int wheel=0;wheel<5;wheel++){
     for(int station=0;station<4;station++){
       for(int sector=0;sector<14;sector++){
	 TString hName = Form("Wh%d_St%d_Sec%d",wheel-2,station+1,sector+1);
	 TString hNamef = Form("Wh%d_St%d_Sec%d",wheel-2,station+1,sector+1);
	 //cout << hName << endl;
	 HL[wheel][station][sector] = new TH1F(Form("%s_HL",hName.Data()),hName,120,-80,40); 
	 HH[wheel][station][sector] = new TH1F(Form("%s_HH",hName.Data()),hName,120,-80,40); 
	 div[wheel][station][sector]= new TH1F(Form("%s_div",hName.Data()),hName,120,-80,40);   
	 HL_fold[wheel][station][sector] = new TH1F(Form("%s_HL_fold",hNamef.Data()),hName,50,0,25); 
	 HH_fold[wheel][station][sector] = new TH1F(Form("%s_HH_fold",hNamef.Data()),hName,50,0,25); 
	 div_fold[wheel][station][sector]= new TH1F(Form("%s_div_fold",hNamef.Data()),hName,50,0,25); 
       }
     }
   }

   
}

Bool_t histograms_HH_and_HL::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef histograms_HH_and_HL_cxx
