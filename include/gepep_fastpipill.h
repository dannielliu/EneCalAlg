//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 23 09:54:31 2014 by ROOT version 5.34/19
// from TTree gepep_fastpipill/ks N-Tuple example
// found on file: data_Rvalue_pipill_1.root
//////////////////////////////////////////////////////////

#ifndef gepep_fastpipill_h
#define gepep_fastpipill_h

#include "TROOT.h"
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
struct Psip;

class gepep_fastpipill {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           rec;
   Int_t           evttag;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Double_t        mcgpx[3];
   Double_t        mcgpy[3];
   Double_t        mcgpz[3];
   Double_t        mcge[3];
   Double_t        mcgid[3];
   Int_t           nGammatch;
   Int_t           ncharg;
   Int_t           ntot;
   Int_t           nneu;
   Int_t           ngch;
   Int_t           ngam;
   Int_t           npi0;
   Int_t           netap2;
   Double_t        delang[30];   //[nneu]
   Double_t        delphi[30];   //[nneu]
   Double_t        delthe[30];   //[nneu]
   Int_t           npart[30];   //[nneu]
   Int_t           nemchits[30];   //[nneu]
   Double_t        module[30];   //[nneu]
   Double_t        x[30];   //[nneu]
   Double_t        y[30];   //[nneu]
   Double_t        z[30];   //[nneu]
   Double_t        px[30];   //[nneu]
   Double_t        py[30];   //[nneu]
   Double_t        pz[30];   //[nneu]
   Double_t        theta[30];   //[nneu]
   Double_t        phi[30];   //[nneu]
   Double_t        dx[30];   //[nneu]
   Double_t        dy[30];   //[nneu]
   Double_t        dz[30];   //[nneu]
   Double_t        dtheta[30];   //[nneu]
   Double_t        dphi[30];   //[nneu]
   Double_t        energy[30];   //[nneu]
   Double_t        dE[30];   //[nneu]
   Double_t        eSeed[30];   //[nneu]
   Double_t        nSeed[30];   //[nneu]
   Double_t        e3x3[30];   //[nneu]
   Double_t        e5x5[30];   //[nneu]
   Double_t        secondMoment[30];   //[nneu]
   Double_t        latMoment[30];   //[nneu]
   Double_t        a20Moment[30];   //[nneu]
   Double_t        a42Moment[30];   //[nneu]
   Double_t        getTime[30];   //[nneu]
   Double_t        getEAll[30];   //[nneu]
   Double_t        mpi0[80];   //[npi0]
   Double_t        chisq1cpi0[80];   //[npi0]
   Int_t           ig1pi0[80];   //[npi0]
   Int_t           ig2pi0[80];   //[npi0]
   Double_t        chisq6c;
   Double_t        chisq4c;
   Double_t        chisq3c;
   Double_t        gpx4[3];
   Double_t        gpy4[3];
   Double_t        gpz4[3];
   Double_t        ge4[3];
   Double_t        lepx4[2];// lepton 4-momentum px [0] for + [1] for -
   Double_t        lepy4[2];
   Double_t        lepz4[2];
   Double_t        lee4[2];
   Double_t        pipx4[2];
   Double_t        pipy4[2];
   Double_t        pipz4[2];
   Double_t        pie4[2];
   Double_t        ggm4;
   Double_t        kkm4;
   Double_t        kpi01m4;
   Double_t        kpi02m4;
   Double_t        kkpi0m4;
   Double_t        kkpipim4;
   Double_t        kkpipipi0m4;
   Double_t        Recoilmass;
   Double_t        llm4;
   Double_t        pipim4;
   Double_t        llpipi1m4;
   Double_t        llpipi2m4;
   Double_t        llpipi3m4;
   Double_t        llpipi4m4;
   Double_t        gpipim4;
   Double_t        decay_ee;
   Double_t        hepp;
   Double_t        hepm;
   Double_t        emcp;
   Double_t        emcm;
   Double_t        mucp;
   Double_t        mucm;
   Double_t        zcpm4;
   Double_t        zcmm4;
   Double_t        ecms;
   Double_t        angle1;
   Double_t        angle2;
   Double_t        angle3;
   Double_t        angle4;
   Double_t        psipm4;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_rec;   //!
   TBranch        *b_evttag;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_mcgpx;   //!
   TBranch        *b_mcgpy;   //!
   TBranch        *b_mcgpz;   //!
   TBranch        *b_mcge;   //!
   TBranch        *b_mcgid;   //!
   TBranch        *b_nGammatch;   //!
   TBranch        *b_ncharg;   //!
   TBranch        *b_ntot;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_ngch;   //!
   TBranch        *b_ngam;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_netap2;   //!
   TBranch        *b_delang;   //!
   TBranch        *b_delphi;   //!
   TBranch        *b_delthe;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_nemchits;   //!
   TBranch        *b_module;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_dx;   //!
   TBranch        *b_dy;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_dtheta;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_dE;   //!
   TBranch        *b_eSeed;   //!
   TBranch        *b_nSeed;   //!
   TBranch        *b_e3x3;   //!
   TBranch        *b_e5x5;   //!
   TBranch        *b_secondMoment;   //!
   TBranch        *b_latMoment;   //!
   TBranch        *b_a20Moment;   //!
   TBranch        *b_a42Moment;   //!
   TBranch        *b_getTime;   //!
   TBranch        *b_getEAll;   //!
   TBranch        *b_mpi0;   //!
   TBranch        *b_chisq1cpi0;   //!
   TBranch        *b_ig1pi0;   //!
   TBranch        *b_ig2pi0;   //!
   TBranch        *b_chisq6c;   //!
   TBranch        *b_chisq4c;   //!
   TBranch        *b_chisq3c;   //!
   TBranch        *b_gpx4;   //!
   TBranch        *b_gpy4;   //!
   TBranch        *b_gpz4;   //!
   TBranch        *b_ge4;   //!
   TBranch        *b_lepx4;   //!
   TBranch        *b_lepy4;   //!
   TBranch        *b_lepz4;   //!
   TBranch        *b_lee4;   //!
   TBranch        *b_pipx4;   //!
   TBranch        *b_pipy4;   //!
   TBranch        *b_pipz4;   //!
   TBranch        *b_pie4;   //!
   TBranch        *b_ggm4;   //!
   TBranch        *b_kkm4;   //!
   TBranch        *b_kpi01m4;   //!
   TBranch        *b_kpi02m4;   //!
   TBranch        *b_kkpi0m4;   //!
   TBranch        *b_kkpipim4;   //!
   TBranch        *b_kkpipipi0m4;   //!
   TBranch        *b_Recoilmass;   //!
   TBranch        *b_llm4;   //!
   TBranch        *b_pipim4;   //!
   TBranch        *b_llpipi1m4;   //!
   TBranch        *b_llpipi2m4;   //!
   TBranch        *b_llpipi3m4;   //!
   TBranch        *b_llpipi4m4;   //!
   TBranch        *b_gpipim4;   //!
   TBranch        *b_decay_ee;   //!
   TBranch        *b_hepp;   //!
   TBranch        *b_hepm;   //!
   TBranch        *b_emcp;   //!
   TBranch        *b_emcm;   //!
   TBranch        *b_mucp;   //!
   TBranch        *b_mucm;   //!
   TBranch        *b_zcpm4;   //!
   TBranch        *b_zcmm4;   //!
   TBranch        *b_ecms;   //!
   TBranch        *b_angle1;   //!
   TBranch        *b_angle2;   //!
   TBranch        *b_angle3;   //!
   TBranch        *b_angle4;   //!
   TBranch        *b_psipm4;   //!

   gepep_fastpipill(TTree *tree=0);
   virtual ~gepep_fastpipill();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   bool     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void FitSpe(std::vector<Psip> &evts, const char* namesfx);
};

#endif

