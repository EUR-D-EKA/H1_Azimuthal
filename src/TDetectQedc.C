#include <iostream>
#include <iomanip>
#include "TDetectQedc.h"
#include "H1Mods/H1PartMCArrayPtr.h"
#include "H1Mods/H1PartMC.h"

using namespace std;

TDetectQedc::TDetectQedc(H1PartMCArrayPtr &mcpart) {
   fIelectron=-2;
   fIelectronNoFSR=-2;
   fNisr=0;
   fNfsr=0;
   fIphotonQedc=-2;
   fELIE=kELIE_DEFAULT;
   fAGMA=kAGMA_DEFAULT;
   fAGMI=kAGMI_DEFAULT;
   fELIG=kELIG_DEFAULT;
   fPTMA=kPTMA_DEFAULT;
   fVMIN=kVMIN_DEFAULT;
   fWMIN=kWMIN_DEFAULT;
   fWMAX=kWMAX_DEFAULT;
   fACO =kACO_DEFAULT;
   fNcluster=0;
   fNstring=0;
   fMstring=0.0;
   fMcluster=0.0;
   fMesonPdg=0;
   Int_t vmParent=-1;
   fElectron.SetXYZT(0.,0.,0.,0.);
   for(Int_t i=0;i<PHOTON_NTYPE;i++) {
      fPhoton[i].SetXYZT(0.,0.,0.,0.);
      fCutsQedc[i]=0;
   }
   fEbeam.SetXYZT(0.,0.,0.,0.);
   fEPbeam.SetXYZT(0.,0.,0.,0.);
   fW.SetXYZT(0.,0.,0.,0.);
   fInvis.SetXYZT(0.,0.,0.,0.);

   // special code for DVCS
   Int_t nDVCS=0;
   for(Int_t i=0;(i<mcpart.GetEntries())&&(i<=6);i++) {
      int pdg=abs(mcpart[i]->GetPDG());
      if(mcpart[i]->GetStatus()!=3) continue;
      if((i==0)&&(pdg!=2212)) continue;
      if((i==1)&&(pdg!=11)) continue;
      if((i==2)&&(pdg!=11)) continue;
      if((i==3)&&(pdg!=22)) continue;
      if((i==4)&&(pdg!=11)) continue;
      if((i==6)&&(pdg!=22)) continue;
      nDVCS++;
   }
   if(nDVCS==7) {
      // DVCS Monte Carlo
      fEPbeam = mcpart[0]->GetFourVector()+mcpart[1]->GetFourVector();
      fEbeam = mcpart[1]->GetFourVector();
      fPhoton[PHOTON_ISR] = mcpart[3]->GetFourVector();
      fNisr++;
      fElectron = mcpart[4]->GetFourVector();
      fElectronNoFSR = mcpart[4]->GetFourVector();
      fIelectronNoFSR = 4;
      fIelectron = 4;
      fPhoton[PHOTON_QEDC]=mcpart[6]->GetFourVector();
      fIphotonQedc=6;
   } else {
      // non-DVCS Monte Carlo
   for(Int_t i=0;i<mcpart.GetEntries();i++) {
      H1PartMC *p=mcpart[i];
      int pdg=abs(p->GetPDG());
      int status=p->GetStatus();
      if(status<0) continue;
      int parent_status=-2;
      int parent_pdg=-1;
      if(p->GetMother1()>=0) {
         parent_status=mcpart[p->GetMother1()]->GetStatus();
         parent_pdg=abs(mcpart[p->GetMother1()]->GetPDG());
      }
      if(status==201) {
         fEPbeam += p->GetFourVector();
         if(pdg==11) fEbeam= p->GetFourVector();
      }
      if(pdg==11) {
         // electron or positron
         if((status==0) && (fIelectron<0)) {
            if((parent_status==201)|| // DJANGO and COMPTON without ISR
               (p->GetMother1()==fIelectronNoFSR) // COMPTON FSR
               ) {
               fElectron=p->GetFourVector();
               if(fIelectronNoFSR<0) {
                  // no FSR, store electron here
                  fElectronNoFSR=p->GetFourVector();
               }
               fIelectron=i;
            }
         } else if((status==2) && (fIelectronNoFSR<0) &&
                   (parent_status==201)) {
            fElectronNoFSR=p->GetFourVector();
            fIelectronNoFSR=i;
         }
      } else if(pdg==22) {
         // photon
         if(((status==202)&&(fIelectron==p->GetMother1())) // DJANGO
            || 
            ((status==0)&&(fIelectronNoFSR==p->GetMother1()))) { // COMPTON FSR
            fPhoton[PHOTON_FSR] += p->GetFourVector();
            fNfsr++;
         } else if((status==202)&&(parent_status==201)) {
            fPhoton[PHOTON_ISR] += p->GetFourVector();
            fNisr++;
         } else if((status==202)||
                   ((status==0) && (parent_pdg<20))) {
            if(fIphotonQedc<0) {
               fPhoton[PHOTON_QEDC]=p->GetFourVector();
               fIphotonQedc=i;
            }
         }
      }
      // count number of clusters and strings
      if(pdg==91) {
         fNcluster++;
         fMcluster += p->GetMass();
      }
      if(pdg==92) {
         fNstring++;
         fMstring += p->GetMass();
      }
      // detect DIFFVM meson: photon going into a single VM
      if((status==3)&&(pdg>100)&&
         (parent_status==3)&&(parent_pdg==22)&&
         (p->GetMother2()==-1)) {
         // check whether there are single daughters only
         if((mcpart[p->GetMother1()]->GetDaughter2()==i)&&
            (p->GetDaughter1()==p->GetDaughter2())) {
            // this is the VM from the documentation line
            vmParent=i;
         }
      }
      if((status==2)&&(pdg>100)&&(vmParent>=0)&&(vmParent==p->GetMother1())&&
         (pdg==parent_pdg)) {
         fMesonPdg=p->GetPDG();
         fMeson=p->GetFourVector();
      }
      // detect DJANGO meson here
      if((!fMesonPdg)&&
         (pdg==91)&&
         (p->GetMother1()>=0)&&(p->GetMother2()>=0)) {
         H1PartMC *m1=mcpart[p->GetMother1()];
         H1PartMC *m2=mcpart[p->GetMother2()];
         if((m1->GetPDG()+m2->GetPDG()==0)&&
            (m1->GetDaughter2()==i)&&
            (m2->GetDaughter2()==i)) {
            // cluster made from q-qbar pair
            if((p->GetDaughter1()==p->GetDaughter2())&&
               (p->GetDaughter1()>=0)) {
               H1PartMC *vm=mcpart[p->GetDaughter1()];
               fMesonPdg=vm->GetPDG();
               fMeson=vm->GetFourVector();
            }
         }
      }
   }
   }
   // check cuts of all pairings
   for(Int_t i=0;i<PHOTON_NTYPE;i++) {
      fCutsQedc[i]=GetQedcCuts(fElectron,fPhoton[i]);
   }
   // calculate more 4-vectors
   fW=fEPbeam-fElectron -fPhoton[PHOTON_QEDC] -fPhoton[PHOTON_FSR];

   // mass of all final-state particles
   // which are neither the electron nor the QEDC photon
   fInvis=fW;
   // mass of the outgoing hadronic system
   fW -= fPhoton[PHOTON_ISR];
}

void TDetectQedc::Print(H1PartMCArrayPtr &mcpart) {
   cout<<"TDetectQedc::Print npart="<<mcpart.GetEntries()
       <<" nString="<<fNstring
       <<" nCluster="<<fNcluster
       <<"\n";
   if(fMesonPdg) {
      cout<<"Meson pdg="<<fMesonPdg
          <<" Pt="<<fMeson.Pt()
          <<" Pz="<<fMeson.Pz()
          <<" E="<<fMeson.E()
          <<" M="<<fMeson.M()<<"\n";
   }
   double tfsr=(fEbeam-fPhoton[PHOTON_ISR]-fPhoton[PHOTON_QEDC]-fPhoton[PHOTON_FSR]-fElectron).M();
   double tnofsr=(fEbeam-fPhoton[PHOTON_ISR]-fPhoton[PHOTON_QEDC]-fElectronNoFSR).M();
   cout<<"t after FSR: sqrt(t)="<<tfsr <<"\n";
   cout<<"t before FSR: sqrt(t)="<<tnofsr <<"\n";
   for(int i=0;i<mcpart.GetEntries();i++) {
      H1PartMC *p=mcpart[i];
      int pdg=abs(p->GetPDG());
      int status=p->GetStatus();
      //if(status<0) continue;
      int parent_status=-2;
      if(p->GetMother1()>=0) {
         parent_status=mcpart[p->GetMother1()]->GetStatus();
      }
      /* cout<<setw(3)<<i<<" "<<setw(5)<<pdg<<" "<<setw(3)<<status
         <<" "<<setw(3)<<p->GetMother1()<<" "<<setw(3)<<p->GetMother2()
         <<" "<<setw(3)<<p->GetDaughter1()<<" "<<setw(3)<<p->GetDaughter2()
         <<" ("<<p->GetE()<<","<<p->GetPt()<<","<<p->GetPz()<<")"
         <<"\n"; */
      if(i==fIelectronNoFSR) cout<<" --- Next is eqedcNoFSR ---\n";
      if(i==fIelectron) cout<<" --- Next is eqedc ---\n";
      if(i==fIphotonQedc) cout<<" --- Next is gqedc ---\n";
      if(pdg==22) {
         if(((status==202)&&(fIelectron==p->GetMother1()))||
            ((status==0)&&(fIelectronNoFSR==p->GetMother1()))) {
            cout<<" --- Next is FSR ---\n";
         } else if((status==202)&&(parent_status==201)) {
            cout<<" --- Next is ISR ---\n";
         }
      }
      if(!i) {
         p->Print("HEAD");
      }
      cout<<setw(3)<<i<<setw(6)<<p->GetPDG()<<" ";
      p->Print();
   }
}

Bool_t TDetectQedc::IsElectronFound(void) const {
   return fIelectron>=0;
}

Bool_t TDetectQedc::IsPhotonQedcFound(void) const {
   return fIphotonQedc>=0;
}

Float_t const TDetectQedc::kELIE_DEFAULT=1.5;
Float_t const TDetectQedc::kAGMA_DEFAULT=178.;
Float_t const TDetectQedc::kAGMI_DEFAULT=3.6;
Float_t const TDetectQedc::kELIG_DEFAULT=1.5;
Float_t const TDetectQedc::kPTMA_DEFAULT=20.;
Float_t const TDetectQedc::kVMIN_DEFAULT=15.;
Float_t const TDetectQedc::kWMIN_DEFAULT=1.5;
Float_t const TDetectQedc::kWMAX_DEFAULT=310.;
Float_t const TDetectQedc::kACO_DEFAULT=50.;

Int_t TDetectQedc::GetQedcCuts
(TLorentzVector const &electron,TLorentzVector const &photon) const {
   Int_t r=0;
   if(electron.E()<fELIE)  r |= 0x00000001;
   Double_t theta=electron.Theta()*180./M_PI;
   if(theta>fAGMA)         r |= 0x00000002;
   if(theta<fAGMI)         r |= 0x00000004;
   if(photon.E()<fELIG)    r |= 0x00000010;
   theta=photon.Theta()*180./M_PI;
   if(theta>fAGMA)         r |= 0x00000020;
   if(theta<fAGMI)         r |= 0x00000040;
   TLorentzVector s=electron+photon;
   if(s.Pt()>fPTMA)        r |= 0x00000080;
   if(s.E()<fVMIN)         r |= 0x00000100;
   if(s.M()<fWMIN)         r |= 0x00000200;
   if(s.M()>fWMAX)         r |= 0x00000400;
   Double_t dphi=fabs(remainder(electron.Phi()+M_PI-photon.Phi(),2.*M_PI)
                      *180./M_PI);
   if(dphi>fACO)           r |= 0x00000800;
   return r;
}

Bool_t TDetectQedc::IsQedcEvent(void) const {
   for(Int_t i=0;i<PHOTON_NTYPE;i++) {
      // any of the pairs inside all QEDC cuts?
      if(!fCutsQedc[i]) return kTRUE;
   }
   // none of the pairs matches all QEDC cuts
   return kFALSE;
}
