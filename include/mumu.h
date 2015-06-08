//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May 31 16:05:42 2015 by ROOT version 5.34/28
// from TTree mumu/ks N-Tuple example
// found on file: Rvalue_mumu_3850.root
//////////////////////////////////////////////////////////

#ifndef mumu_h
#define mumu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class mumu {
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
   Double_t        lepx4[2];
   Double_t        lepy4[2];
   Double_t        lepz4[2];
   Double_t        lee4[2];
   Double_t        pipx4[2];
   Double_t        pipy4[2];
   Double_t        pipz4[2];
   Double_t        pie4[2];
   Double_t        llm4;
   Double_t        pipim4;
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
   TBranch        *b_lepx4;   //!
   TBranch        *b_lepy4;   //!
   TBranch        *b_lepz4;   //!
   TBranch        *b_lee4;   //!
   TBranch        *b_pipx4;   //!
   TBranch        *b_pipy4;   //!
   TBranch        *b_pipz4;   //!
   TBranch        *b_pie4;   //!
   TBranch        *b_llm4;   //!
   TBranch        *b_pipim4;   //!
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

   mumu(TTree *tree=0);
   virtual ~mumu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mumu_cxx
mumu::mumu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Rvalue_mumu_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Rvalue_mumu_3850.root");
      }
      f->GetObject("mumu",tree);

   }
   Init(tree);
}

mumu::~mumu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mumu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mumu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mumu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("rec", &rec, &b_rec);
   fChain->SetBranchAddress("evttag", &evttag, &b_evttag);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("mcgpx", mcgpx, &b_mcgpx);
   fChain->SetBranchAddress("mcgpy", mcgpy, &b_mcgpy);
   fChain->SetBranchAddress("mcgpz", mcgpz, &b_mcgpz);
   fChain->SetBranchAddress("mcge", mcge, &b_mcge);
   fChain->SetBranchAddress("mcgid", mcgid, &b_mcgid);
   fChain->SetBranchAddress("nGammatch", &nGammatch, &b_nGammatch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("ntot", &ntot, &b_ntot);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ngam", &ngam, &b_ngam);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("netap2", &netap2, &b_netap2);
   fChain->SetBranchAddress("delang", delang, &b_delang);
   fChain->SetBranchAddress("delphi", delphi, &b_delphi);
   fChain->SetBranchAddress("delthe", delthe, &b_delthe);
   fChain->SetBranchAddress("npart", npart, &b_npart);
   fChain->SetBranchAddress("nemchits", nemchits, &b_nemchits);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dx", dx, &b_dx);
   fChain->SetBranchAddress("dy", dy, &b_dy);
   fChain->SetBranchAddress("dz", dz, &b_dz);
   fChain->SetBranchAddress("dtheta", dtheta, &b_dtheta);
   fChain->SetBranchAddress("dphi", dphi, &b_dphi);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("dE", dE, &b_dE);
   fChain->SetBranchAddress("eSeed", eSeed, &b_eSeed);
   fChain->SetBranchAddress("nSeed", nSeed, &b_nSeed);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("secondMoment", secondMoment, &b_secondMoment);
   fChain->SetBranchAddress("latMoment", latMoment, &b_latMoment);
   fChain->SetBranchAddress("a20Moment", a20Moment, &b_a20Moment);
   fChain->SetBranchAddress("a42Moment", a42Moment, &b_a42Moment);
   fChain->SetBranchAddress("getTime", getTime, &b_getTime);
   fChain->SetBranchAddress("getEAll", getEAll, &b_getEAll);
   fChain->SetBranchAddress("mpi0", mpi0, &b_mpi0);
   fChain->SetBranchAddress("chisq1cpi0", chisq1cpi0, &b_chisq1cpi0);
   fChain->SetBranchAddress("ig1pi0", ig1pi0, &b_ig1pi0);
   fChain->SetBranchAddress("ig2pi0", ig2pi0, &b_ig2pi0);
   fChain->SetBranchAddress("lepx4", lepx4, &b_lepx4);
   fChain->SetBranchAddress("lepy4", lepy4, &b_lepy4);
   fChain->SetBranchAddress("lepz4", lepz4, &b_lepz4);
   fChain->SetBranchAddress("lee4", lee4, &b_lee4);
   fChain->SetBranchAddress("pipx4", pipx4, &b_pipx4);
   fChain->SetBranchAddress("pipy4", pipy4, &b_pipy4);
   fChain->SetBranchAddress("pipz4", pipz4, &b_pipz4);
   fChain->SetBranchAddress("pie4", pie4, &b_pie4);
   fChain->SetBranchAddress("llm4", &llm4, &b_llm4);
   fChain->SetBranchAddress("pipim4", &pipim4, &b_pipim4);
   fChain->SetBranchAddress("decay_ee", &decay_ee, &b_decay_ee);
   fChain->SetBranchAddress("hepp", &hepp, &b_hepp);
   fChain->SetBranchAddress("hepm", &hepm, &b_hepm);
   fChain->SetBranchAddress("emcp", &emcp, &b_emcp);
   fChain->SetBranchAddress("emcm", &emcm, &b_emcm);
   fChain->SetBranchAddress("mucp", &mucp, &b_mucp);
   fChain->SetBranchAddress("mucm", &mucm, &b_mucm);
   fChain->SetBranchAddress("zcpm4", &zcpm4, &b_zcpm4);
   fChain->SetBranchAddress("zcmm4", &zcmm4, &b_zcmm4);
   fChain->SetBranchAddress("ecms", &ecms, &b_ecms);
   fChain->SetBranchAddress("angle1", &angle1, &b_angle1);
   fChain->SetBranchAddress("angle2", &angle2, &b_angle2);
   fChain->SetBranchAddress("angle3", &angle3, &b_angle3);
   fChain->SetBranchAddress("angle4", &angle4, &b_angle4);
   fChain->SetBranchAddress("psipm4", &psipm4, &b_psipm4);
   Notify();
}

Bool_t mumu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mumu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mumu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mumu_cxx
