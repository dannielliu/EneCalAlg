#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "Ks0Alg_sim.h"

int main(int argc, char **argv)
{
  char filename[1000];
  char outdir[1000] = {"/Volumes/data2/newFiles"};
  if (argc==1){
    std::cout<<"Error: you should add root file name as input arguements"<<std::endl;
	return -1;
  }
  TFile *file;
  int fileNo = argc-1;

    sprintf(filename,"%s",argv[1]);
    file=new TFile(filename);
    if(!file){
      std::cout<<"can not open file "<<filename<<std::endl;
      return false;
    }
  
  
  char treename[100] = {"Ks_info"};
  TTree *tree = new TTree();
  file->GetObject(treename,tree);
  if(!tree){
    std::cout<<"can not find tree "<<treename<<std::endl;
    return false;
  }
  char newname[1000];
  int idx=0;
  for (int i=0; i<1000;i++){
    if (filename[i] == '/') idx=i;
	if (filename[i] == '\0') break;
  }
  std::cout<<&filename[idx+1]<<std::endl;
  sprintf(newname,"%s/new_%s",outdir,&filename[idx+1]);
  std::cout<<newname<<std::endl;
  TFile *newfile = new TFile(newname,"recreate");
  TTree *newtree = new TTree("newtree",treename);
  
  Ks0Alg alg;
  alg.Init(tree);
  alg.Link(newtree);
  newtree->SetAutoSave(300000);
  alg.Loop(newtree);

  newtree->Write();
  newfile->Close();
  file->Close();
  delete newtree;
  delete tree;
  delete newfile;
  delete file;

  return 0;
}

#define Ks0Alg_cxx
#include <TH2.h>
#include <vector>

void Ks0Alg::Loop(TTree *tree)
{
//   In a ROOT session, you can do:
//      Root > .L Ks0Alg.C
//      Root > Ks0Alg t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

  
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      tree->Fill();
	  if (jentry%10000 == 0)
		std::cout<<jentry<<" events......"<<std::endl;
   }// loop data end 

   return;
}



#ifdef Ks0Alg_cxx
Ks0Alg::Ks0Alg(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ksto2pi_583_140226_0035862.root");
    //if (!f) {
    //   f = new TFile("Ksto2pi_583_140226_0035862.root");
    //}
    //tree = (TTree*)gDirectory->Get("Ks_info");
      std::cout<<"Waring: no input tree when construct"<<std::endl;
	  return;
   }
   Init(tree);
}

Ks0Alg::~Ks0Alg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ks0Alg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ks0Alg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ks0Alg::Init(TTree *tree)
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

   fChain->SetBranchAddress("rec_truth_Mks", &rec_truth_Mks, &b_rec_truth_Mks);
   fChain->SetBranchAddress("rec_truth_Pks", &rec_truth_Pks, &b_rec_truth_Pks);
   fChain->SetBranchAddress("rec_truth_Eks", &rec_truth_Eks, &b_rec_truth_Eks);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("nKsCan", &nKsCan, &b_nKsCan);
   fChain->SetBranchAddress("m_chisq0", m_chisq0, &b_m_chisq0);
   fChain->SetBranchAddress("m_chisqvtxnd", m_chisqvtxnd, &b_m_chisqvtxnd);
   fChain->SetBranchAddress("m_decayL", m_decayL, &b_m_decayL);
   fChain->SetBranchAddress("m_decayLerr", m_decayLerr, &b_m_decayLerr);
   fChain->SetBranchAddress("m_ctau", m_ctau, &b_m_ctau);
   fChain->SetBranchAddress("Ks_id", Ks_id, &b_Ks_id);
   fChain->SetBranchAddress("Mpippim", Mpippim, &b_Mpippim);
   fChain->SetBranchAddress("Ppippim", Ppippim, &b_Ppippim);
   fChain->SetBranchAddress("Epippim", Epippim, &b_Epippim);
   fChain->SetBranchAddress("Thepippim", Thepippim, &b_Thepippim);
   fChain->SetBranchAddress("Phipippim", Phipippim, &b_Phipippim);
   fChain->SetBranchAddress("Ppip", Ppip, &b_Ppip);
   fChain->SetBranchAddress("Ppim", Ppim, &b_Ppim);
   fChain->SetBranchAddress("Thepip", Thepip, &b_Thepip);
   fChain->SetBranchAddress("Thepim", Thepim, &b_Thepim);
   fChain->SetBranchAddress("Phipip", Phipip, &b_Phipip);
   fChain->SetBranchAddress("Phipim", Phipim, &b_Phipim);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   Notify();
}

void Ks0Alg::Link(TTree *tree)
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
   //fChain = tree;
   //fCurrent = -1;
   //fChain->SetMakeClass(1);

   tree->Branch("rec_truth_Mks", &rec_truth_Mks, "rec_truth_Mks/D");
   tree->Branch("rec_truth_Pks", &rec_truth_Pks, "rec_truth_Pks/D");
   tree->Branch("rec_truth_Eks", &rec_truth_Eks, "rec_truth_Eks/D");
   tree->Branch("indexmc", &indexmc, "indexmc/I");
   tree->Branch("pdgid", pdgid, "pdgid[100]/I");
   tree->Branch("motheridx", motheridx, "motheridx[100]/I");
   tree->Branch("npip", &npip, "npip/I");
   tree->Branch("npim", &npim, "npim/I");
   tree->Branch("nKsCan", &nKsCan, "nKsCan/I");
   tree->Branch("m_chisq0", m_chisq0, "m_chisq0[100]/D");
   tree->Branch("m_chisqvtxnd", m_chisqvtxnd, "m_chisqvtxnd[100]/D");
   tree->Branch("m_decayL", m_decayL, "m_decayL[100]/D");
   tree->Branch("m_decayLerr", m_decayLerr, "m_decayLerr[100]/D");
   tree->Branch("m_ctau", m_ctau, "m_ctau[100]/D");
   tree->Branch("Ks_id", Ks_id, "Ks_id[100]/I");
   tree->Branch("Mpippim", Mpippim, "Mpippim[100]/D");
   tree->Branch("Ppippim", Ppippim, "Ppippim[100]/D");
   tree->Branch("Epippim", Epippim, "Epippim[100]/D");
   tree->Branch("Thepippim", Thepippim, "Thepippim[100]/D");
   tree->Branch("Phipippim", Phipippim, "Phipippim[100]/D");
   tree->Branch("Ppip", Ppip, "Ppip[100]/D");
   tree->Branch("Ppim", Ppim, "Ppim[100]/D");
   tree->Branch("Thepip", Thepip, "Thepip[100]/D");
   tree->Branch("Thepim", Thepim, "Thepim[100]/D");
   tree->Branch("Phipip", Phipip, "Phipip[100]/D");
   tree->Branch("Phipim", Phipim, "Phipim[100]/D");
   tree->Branch("pippx", pippx, "pippx[100]/D");
   tree->Branch("pippy", pippy, "pippy[100]/D");
   tree->Branch("pippz", pippz, "pippz[100]/D");
   tree->Branch("pipe", pipe, "pipe[100]/D");
   tree->Branch("pimpx", pimpx, "pimpx[100]/D");
   tree->Branch("pimpy", pimpy, "pimpy[100]/D");
   tree->Branch("pimpz", pimpz, "pimpz[100]/D");
   tree->Branch("pime", pime, "pime[100]/D");
}


Bool_t Ks0Alg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ks0Alg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ks0Alg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ks0Alg_cxx
