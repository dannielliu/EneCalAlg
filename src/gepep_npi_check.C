#define gepep_npi_cxx
#include "gepep_npi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
extern std::string outputdir;

void gepep_npi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_npi.C
//      Root > gepep_npi t
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

   char fname[1000];
   double mpi = 0.13957;
   double mD0 = 1.86484;
   double mass;
   double costheta1,costheta2;
   sprintf(fname,"%s/plot_npi.root",outputdir.c_str());
   TFile *f = new TFile(fname,"recreate");
   TTree *vars = new TTree("vars","vars");
   vars->Branch("mass",&mass,"mass/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   NPi evt;

   Long64_t nbytes = 0, nb = 0;
   nentries = nentries > 100000? 100000: nentries;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  //for (int ip=0;ip<npip;ip++)
	  //for (int im=0;im<npim;im++){
	  //  evt.SetVal(pippx[ip],pippy[ip],pippz[ip],pimpx[im],pimpy[im],pimpz[im]);
	  //  mass = evt.InvMass();
	  //  costheta1 = evt.GetCostheta1();
	  //  costheta2 = evt.GetCostheta2();
	  //  vars->Fill();
	  //
	  //}
	//if (npip==3 && npim==3){
	//  evt.Clear();
	//  evt.AddParticle(pippx[0],pippy[0],pippz[0]);
	//  evt.AddParticle(pippx[1],pippy[1],pippz[1]);
	//  evt.AddParticle(pippx[2],pippy[2],pippz[2]);
	//  evt.AddParticle(pimpx[0],pimpy[0],pimpz[0]);
	//  evt.AddParticle(pimpx[1],pimpy[1],pimpz[1]);
	//  evt.AddParticle(pimpx[2],pimpy[2],pimpz[2]);
	//  mass = evt.InvMass();
	//  vars->Fill();

	//}
	//else{
	//  //std::cout<<"event id "<<jentry<<" is not a 4 pi event"<<std::endl;
	//  continue;
	//}
	  
	  if (npip==0 || npim==0) continue;
	  NPi tmpEvt;
	  double massgap =10;
	  for (int i=0; i<npip; i++){
	  for (int j=0; j<npim; j++){
	    tmpEvt.AddParticle(pippx[i],pippy[i],pippz[i]);
	    tmpEvt.AddParticle(pimpx[i],pimpy[i],pimpz[i]);
		if (tmpEvt.InvMass()-mD0 > massgap) continue;
		massgap = tmpEvt.InvMass() - mD0;
		evt = tmpEvt;
	  }
	  }
	  mass = evt.InvMass();
	  vars->Fill();

   }
   vars->Write();
   f->Close();
   return;
}
