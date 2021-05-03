///////////////////////////////////////////////////////
// Search for J/Psi on MODS
//
// Author     : Ingo Strauch (strauch@mail.desy.de)
// Created    : 19.12.2000
// Redone by  : Gero Flucke (gero.flucke@desy.de)
//          and Bengt Wessling (bengt.wessling@desy.de)
// Last update: $Date: 2010/10/08 13:23:57 $
//          by: $Author: msteder $
//
///////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <set>

// ROOT includes
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

// for testing
#include "H1Skeleton/H1EventFiller.h"

// H1 OO includes
#include "H1Skeleton/H1Tree.h"
#include "H1Steering/H1StdCmdLine.h"
#include "H1Arrays/H1ArrayF.h"
#include "H1Pointers/H1IntPtr.h"
#include "H1Pointers/H1BytePtr.h"
#include "H1Pointers/H1ShortPtr.h"
#include "H1Pointers/H1FloatPtr.h"
#include "H1Mods/H1PartMCArrayPtr.h"
#include "H1Mods/H1PartMC.h"
#include "H1Mods/H1GetPartMCId.h"
#include "H1Mods/H1SelVertexArrayPtr.h"
#include "H1Mods/H1SelVertex.h"
#include "H1Mods/H1PartCandArrayPtr.h"
#include "H1Mods/H1PartCand.h"
#include "H1Mods/H1PartEm.h"
#include "H1Mods/H1PartSelTrack.h"
#include "H1PhysUtils/H1NuclIACor.h"

#include "H1Geom/H1DetectorStatus.h"
#include "H1Geom/H1DBManager.h"

#include "H1Tracks/H1FSTFittedTrack.h"
#include "H1Tracks/H1FSTFittedTrackArrayPtr.h"

//#include "H1Mods/H1GkiInfoArrayPtr.h"
//#include "H1Mods/H1GkiInfo.h"
#include <TLorentzRotation.h>
#include "H1HadronicCalibration/H1HadronicCalibration.h"


using namespace std;

static double const ME=0.0005109989461;
static double const M_CHARGED_PION=0.13957061;

static double const ELEC_ISOLATION_CONE=0.1;

bool floatEqual(double a,double b) {
   double d1=fabs(a-b);
   double d2=fabs(a+b);
   return (d1<=1.E-5*d2);
}

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());
   return boost;
}

void GetKinematics(TLorentzVector const &ebeam,TLorentzVector const &pbeam,
                   TLorentzVector const &escat,
                   Float_t *x,Float_t *y,Float_t *Q2) {
   TLorentzVector q=ebeam-escat;
   *Q2= -q.Mag2();
   double pq=pbeam.Dot(q);
   *y= pq/pbeam.Dot(ebeam);
   *x= *Q2/2.*pq;
}


struct MyEvent {
   // general information
   Int_t run,evno; // run and event number
   Float_t w; // event weight

   // trigger information
   UInt_t l1l2l3ac[4];
   UInt_t l1l2l3rw[4];
   Float_t trigWeightRW;
   Float_t trigWeightAC;

   // background finders
   UInt_t ibg;
   UInt_t ibgfm;
   UInt_t ibgam;
   UInt_t iqn;

   // event quality .. not complete yet
   Int_t vertexType;
   Float_t vertex[3];
   Float_t beamSpot[2];
   Float_t beamTilt[2];
   Float_t eProtonBeam,eElectronBeam;

   // Monte Carlo information
   Float_t eProtonBeamMC;
   Float_t eElectronBeamMC;
   Float_t simvertex[3];
   Float_t xGKI,yGKI,Q2GKI;

   Float_t elecPxMC,elecPyMC,elecPzMC,elecEMC,elecEradMC; // scattered electron
   Float_t xMC,yMC,Q2MC;
   
   enum {
      nMCtrack_MAX=400
   };
   // if there is no MC info, nMCtrack is set to zero
   Int_t nMCtrackAll;
   Int_t nMCtrack;
   Int_t idMC[nMCtrack_MAX];
#ifdef SAVE_LAB_MOMENTA
   Float_t pxMC[nMCtrack_MAX];
   Float_t pyMC[nMCtrack_MAX];
   Float_t pzMC[nMCtrack_MAX];
#endif
   Float_t ptStarMC[nMCtrack_MAX];
   Float_t etaStarMC[nMCtrack_MAX];
   Float_t phiStarMC[nMCtrack_MAX];
   Float_t log10zMC[nMCtrack_MAX];
   Int_t imatchMC[nMCtrack_MAX];
   
   // reconstructed quantities
   Float_t elecPxREC,elecPyREC,elecPzREC,elecEREC,elecEradREC; // scattered electron
   Float_t xREC,yREC,Q2REC;
   Float_t hfsPxREC,hfsPyREC,hfsPzREC,hfsEREC; // hadronic final state
   enum {
      nRECtrack_MAX=200
   };
   // it there is no scattered electron, no tracks are saved either
   //  (nRECtrack=0)
   Int_t nRECtrackAll;
   Int_t nRECtrack;
   Int_t typeChgREC[nRECtrack_MAX];
#ifdef SAVE_LAB_MOMENTA
   Float_t pxREC[nRECtrack_MAX];
   Float_t pyREC[nRECtrack_MAX];
   Float_t pzREC[nRECtrack_MAX];
#endif
   Float_t ptStarREC[nRECtrack_MAX];
   Float_t etaStarREC[nRECtrack_MAX];
   Float_t phiStarREC[nRECtrack_MAX];
   Float_t log10zREC[nRECtrack_MAX];
   Float_t nucliaREC[nRECtrack_MAX];
   Float_t dmatchREC[nRECtrack_MAX];
   Int_t imatchREC[nRECtrack_MAX];

   // FST tracks
   Int_t nRECfstFitted;
   Int_t nRECfstSelected;

   // auxillary information, not to be saved
   H1PartMC const *partMC[nMCtrack_MAX];
   TVector3 momREC[nRECtrack_MAX];
   TMatrix covREC[nRECtrack_MAX];

};

class DummyFiller : public H1EventFiller {
public:
   virtual void   Fill(H1Event* event) { }
};

int main(int argc, char* argv[]) {
   // parse the command line
   H1StdCmdLine opts;
   opts.Parse(&argc, argv);

   // Load mODS/HAT files
   H1Tree::Instance()->Open();            // this statement must be there!

   //cout << H1Tree::Instance()->SelectHat("NumJPsi>0")
   //     << " events selected " << endl;

   TFile file(opts.GetOutput(), "RECREATE");
   TTree *output=new TTree("properties","properties");
   MyEvent myEvent;
   output->Branch("run",&myEvent.run,"run/I");
   output->Branch("evno",&myEvent.evno,"evno/I");
   output->Branch("w",&myEvent.w,"w/F");
   output->Branch("vertexType",&myEvent.vertexType,"vertexType/I");
   output->Branch("vertex",myEvent.vertex,"vertex[3]/F");
   output->Branch("beamSpot",myEvent.beamSpot,"beamSpot[2]/F");
   output->Branch("beamTilt",myEvent.beamTilt,"beamTilt[2]/F");
   output->Branch("eProtonBeam",&myEvent.eProtonBeam,"eProtonBeam/F");
   output->Branch("eElectronBeam",&myEvent.eElectronBeam,"eElectronBeam/F");
   output->Branch("eProtonBeamMC",&myEvent.eProtonBeamMC,"eProtonBeamMC/F");
   output->Branch("eElectronBeamMC",&myEvent.eElectronBeamMC,"eElectronBeamMC/F");

   output->Branch("l1l2l3ac",myEvent.l1l2l3ac,"l1l2l3ac[4]/i");
   output->Branch("l1l2l3rw",myEvent.l1l2l3ac,"l1l2l3rw[4]/i");
   output->Branch("trigWeightAC",&myEvent.trigWeightAC,"trigWeightAC/F");
   output->Branch("trigWeightRW",&myEvent.trigWeightRW,"trigWeightRW/F");

   output->Branch("ibg",&myEvent.ibg,"ibg/I");
   output->Branch("ibgfm",&myEvent.ibgfm,"ibgfm/I");
   output->Branch("ibgam",&myEvent.ibgam,"ibgam/I");
   output->Branch("iqn",&myEvent.iqn,"iqn/I");

   output->Branch("simvertex",myEvent.simvertex,"simvertex[3]/F");
   output->Branch("elecEradMC",&myEvent.elecEradMC,"elecEradMC/F");
   output->Branch("elecPxMC",&myEvent.elecPxMC,"elecPxMC/F");
   output->Branch("elecPyMC",&myEvent.elecPyMC,"elecPyMC/F");
   output->Branch("elecPzMC",&myEvent.elecPzMC,"elecPzMC/F");
   output->Branch("elecEMC",&myEvent.elecEMC,"elecEMC/F");
   output->Branch("xGKI",&myEvent.xGKI,"xGKI/F");
   output->Branch("yGKI",&myEvent.yGKI,"yGKI/F");
   output->Branch("Q2GKI",&myEvent.Q2GKI,"Q2GKI/F");
   output->Branch("xMC",&myEvent.xMC,"xMC/F");
   output->Branch("yMC",&myEvent.yMC,"yMC/F");
   output->Branch("Q2MC",&myEvent.Q2MC,"Q2MC/F");

   output->Branch("nMCtrackAll",&myEvent.nMCtrackAll,"nMCtrackAll/I");
   output->Branch("nMCtrack",&myEvent.nMCtrack,"nMCtrack/I");
   output->Branch("idMC",myEvent.idMC,"idMC[nMCtrack]/I");
#ifdef SAVE_LAB_MOMENTA
   output->Branch("pxMC",myEvent.pxMC,"pxMC[nMCtrack]/F");
   output->Branch("pyMC",myEvent.pyMC,"pyMC[nMCtrack]/F");
   output->Branch("pzMC",myEvent.pzMC,"pzMC[nMCtrack]/F");
#endif
   output->Branch("ptStarMC",myEvent.ptStarMC,"ptStarMC[nMCtrack]/F");
   output->Branch("etaStarMC",myEvent.etaStarMC,"etaStarMC[nMCtrack]/F");
   output->Branch("phiStarMC",myEvent.phiStarMC,"phiStarMC[nMCtrack]/F");
   output->Branch("log10zMC",myEvent.log10zMC,"log10zMC[nMCtrack]/F");
   output->Branch("imatchMC",myEvent.imatchMC,"imatchMC[nMCtrack]/I");

   output->Branch("elecEradREC",&myEvent.elecEradREC,"elecEradREC/F");
   output->Branch("elecPxREC",&myEvent.elecPxREC,"elecPxREC/F");
   output->Branch("elecPyREC",&myEvent.elecPyREC,"elecPyREC/F");
   output->Branch("elecPzREC",&myEvent.elecPzREC,"elecPzREC/F");
   output->Branch("elecEREC",&myEvent.elecEREC,"elecEREC/F");
   output->Branch("xREC",&myEvent.xREC,"xREC/F");
   output->Branch("yREC",&myEvent.yREC,"yREC/F");
   output->Branch("Q2REC",&myEvent.Q2REC,"Q2REC/F");
   output->Branch("hfsPxREC",&myEvent.hfsPxREC,"hfsPxREC/F");
   output->Branch("hfsPyREC",&myEvent.hfsPyREC,"hfsPyREC/F");
   output->Branch("hfsPzREC",&myEvent.hfsPzREC,"hfsPzREC/F");
   output->Branch("hfsEREC",&myEvent.hfsEREC,"hfsEREC/F");

   output->Branch("nRECtrackAll",&myEvent.nRECtrackAll,"nRECtrackAll/I");
   output->Branch("nRECtrack",&myEvent.nRECtrack,"nRECtrack/I");
   output->Branch("typeChgREC",myEvent.typeChgREC,"typeChgREC[nRECtrack]/I");
#ifdef SAVE_LAB_MOMENTA
   output->Branch("pxREC",myEvent.pxREC,"pxREC[nRECtrack]/F");
   output->Branch("pyREC",myEvent.pyREC,"pyREC[nRECtrack]/F");
   output->Branch("pzREC",myEvent.pzREC,"pzREC[nRECtrack]/F");
#endif
   output->Branch("ptStarREC",myEvent.ptStarREC,"ptStarREC[nRECtrack]/F");
   output->Branch("etaStarREC",myEvent.etaStarREC,"etaStarREC[nRECtrack]/F");
   output->Branch("phiStarREC",myEvent.phiStarREC,"phiStarREC[nRECtrack]/F");
   output->Branch("log10zREC",myEvent.log10zREC,"log10zREC[nRECtrack]/F");
   output->Branch("nucliaREC",myEvent.nucliaREC,"nucliazREC[nRECtrack]/F");
   output->Branch("dmatchREC",myEvent.dmatchREC,"dmatchzREC[nRECtrack]/F");
   output->Branch("imatchREC",myEvent.imatchREC,"imatchzREC[nRECtrack]/I");

   output->Branch("nRECfstFitted",&myEvent.nRECfstFitted,"nRECfstFitted/I");
   output->Branch("nRECfstSelected",&myEvent.nRECfstSelected,
                  "nRECfstSelected/I");

   H1ShortPtr runtype("RunType"); // 0=data, 1=MC, 2=CERN test, 3=CERN MC test
   H1FloatPtr beamx0("BeamX0");          // x position of beam spot (at z=0)
   H1FloatPtr beamy0("BeamY0");          // y position of beam spot (at z=0)

   H1FloatPtr eBeamP("EBeamP"); // proton beam energy from DMIS
   H1FloatPtr eBeamE("EBeamE"); // electron beam energy from DMIS
   H1FloatPtr beamtiltx0("BeamTiltX0");      // x slope of beam axis
   H1FloatPtr beamtilty0("BeamTiltY0");      // y slope of beam axis
   H1IntPtr run("RunNumber");
   H1IntPtr evno("EventNumber");

   H1BytePtr l1l2l3ac("Il1l2l3ac");
   H1BytePtr l1l2rw("Il1l2rw");
   H1BytePtr l1l3rw("Il1l3rw"); 

   H1IntPtr   ibg("Ibg");  // standard array of background finders from OM group (bit packed)
   H1IntPtr  ibgfm("Ibgfm");  // extra finders from OM  (bit packed)
   H1IntPtr  ibgam("Ibgam");  // background finders from DUK group (bit packed)
   H1IntPtr  iqn("Iqn");           // The Lar coherent noise flag

   H1FloatPtr Q2Gki("Q2Gki");
   H1FloatPtr xGki("XGki");
   H1FloatPtr yGki("YGki");
   H1FloatPtr weight1("Weight1");
   H1FloatPtr weight2("Weight2");

   H1FloatPtr genEnElec("GenEnElec"); //   Electron energy, combined with photon for FSR
   H1FloatPtr genPhElec("GenPhElec"); //   Electron phi, combined with photon for FSR
   H1FloatPtr genThElec("GenThElec"); //   Electron theta, combined with photon for FSR

   H1ShortPtr ivtyp("Ivtyp"); // vertex type
   H1SelVertexArrayPtr vertex; // all good vertices
   H1PartCandArrayPtr partCand; // all good tracks

   H1PartMCArrayPtr mcpart;

   H1FSTFittedTrackArrayPtr fstFittedTrack;

   Int_t eventCounter = 0;

   H1HadronicCalibration *hadronicCalibration=H1HadronicCalibration::Instance();
   hadronicCalibration->ApplyHadronicCalibration(H1HadronicCalibration::eHighPtJet);
   hadronicCalibration->ApplyHadronicCalibration(kTRUE);

   // Loop over events
   static int print=10;
   while (gH1Tree->Next() && !opts.IsMaxEvent(eventCounter)) {
      ++eventCounter;
      double w=*weight1 * *weight2;
      // if(*Q2Gki<10.) continue;
      if(print || ((eventCounter %10000)==0))  { 
         cout<<eventCounter
             <<" event "<<*run<<" "<<*evno<<" type="<<*runtype<<" weight="<<w<<"\n";
         if(!print) print=1; // print this event
      }
      myEvent.run=*run;
      myEvent.evno=*evno;
      myEvent.w=w;

      if(*runtype==1) {
         // handle MC information
         if(print) {
            // kinematic variables from GKI bank
            cout<<" xGKI="<<*xGki<<" Q2GKI="<<*Q2Gki
                <<" y="<<*yGki
                <<" w="<<w<<"\n";
         }
         myEvent.xGKI = *xGki;
         myEvent.yGKI = *yGki;
         myEvent.Q2GKI = *Q2Gki;

         H1GetPartMCId mcPartId(&*mcpart);
         mcPartId.Fill();

         TLorentzVector ebeam_MC_lab
            (mcpart[mcPartId.GetIdxBeamElectron()]->GetFourVector());
         TLorentzVector pbeam_MC_lab
            (mcpart[mcPartId.GetIdxBeamProton()]->GetFourVector());

         myEvent.eProtonBeamMC=pbeam_MC_lab.E();
         myEvent.eElectronBeamMC=ebeam_MC_lab.E();

         TLorentzVector escat0_MC_lab
            (mcpart[mcPartId.GetIdxScatElectron()]->GetFourVector());
         // add radiative photon(s) in a cone
         TLorentzVector escatPhot_MC_lab(escat0_MC_lab);
         set<int> isElectron;
         isElectron.insert(mcPartId.GetIdxScatElectron());
         for(int i=0;i<mcpart.GetEntries();i++) {
            H1PartMC *part=mcpart[i];
            int status=part->GetStatus();
            if((status==0 )&&(part->GetPDG()==22)) {
               TLorentzVector p(part->GetFourVector());
               if(p.DeltaR(escat0_MC_lab)<ELEC_ISOLATION_CONE) {
                  // this photon counts with the electron
                  isElectron.insert(i);
                  escatPhot_MC_lab += p;
               }
            }
         }
         myEvent.elecEradMC=escatPhot_MC_lab.E()-escat0_MC_lab.E();
         myEvent.elecPxMC=escatPhot_MC_lab.X();
         myEvent.elecPyMC=escatPhot_MC_lab.Y();
         myEvent.elecPzMC=escatPhot_MC_lab.Z();
         myEvent.elecEMC=escatPhot_MC_lab.E();

         if(print) {
            cout<<"MC scattered electron is made of "<<isElectron.size()<<" particle(s)\n";
         }

         GetKinematics(ebeam_MC_lab,pbeam_MC_lab,escatPhot_MC_lab,
                       &myEvent.xMC,&myEvent.yMC,&myEvent.Q2MC);

         TLorentzRotation boost_MC_HCM
            (BoostToHCM(ebeam_MC_lab,pbeam_MC_lab,escatPhot_MC_lab));
         TLorentzVector q_MC_lab(ebeam_MC_lab-escatPhot_MC_lab);

         // final state particles
         bool haveElectron=false;
         myEvent.nMCtrackAll=0;
         myEvent.nMCtrack=0;
         for(int i=0;i<mcpart.GetEntries();i++) {
            // skip particles counted as electron
            if(isElectron.find(i)!=isElectron.end()) continue;

            H1PartMC *part=mcpart[i];
            if(print) {
               //part->Print();
            }
            int status=part->GetStatus();
            if(status==0) {
               // generator "stable" particles
               if((!haveElectron)&&
                  ((part->GetPDG()==11)||(part->GetPDG()== -11))) {
                  haveElectron=true;               
               } else if(part->GetCharge()!=0.) {
                  // other charged particles
                  TLorentzVector h=part->GetFourVector();
                  double log10z=TMath::Log10((h*pbeam_MC_lab)/(q_MC_lab*pbeam_MC_lab));
                  // boost to hadronic-centre-of-mass frame
                  TLorentzVector hStar = boost_MC_HCM*h;
                  double etaStar=hStar.Eta();
                  double ptStar=hStar.Pt();
                  double phiStar=hStar.Phi();
                  if(print) {
                     cout<<"MCpart "<<myEvent.nMCtrackAll
                         <<" "<<part->GetPDG()
                         <<" etaLab="<<h.Eta()
                         <<" ptLab="<<h.Pt()
                         <<" ptStar="<<ptStar
                         <<" etaStar="<<etaStar
                         <<" phiStar="<<phiStar
                         <<" log10(z)="<<log10z<<"\n";
                  }
                  myEvent.nMCtrackAll++;
                  if(myEvent.nMCtrack<MyEvent::nMCtrack_MAX) {
                     int k=myEvent.nMCtrack;
                     myEvent.idMC[k]=part->GetPDG();
#ifdef SAVE_LAB_MOMENTA
                     myEvent.pxMC[k]=h.X();
                     myEvent.pyMC[k]=h.Y();
                     myEvent.pzMC[k]=h.Z();
#endif
                     myEvent.ptStarMC[k]=hStar.Pt();
                     myEvent.etaStarMC[k]=hStar.Eta();
                     myEvent.phiStarMC[k]=hStar.Phi();
                     myEvent.log10zMC[k]=log10z;
                     myEvent.imatchMC[k]=-1;
                     myEvent.partMC[k]=part;
                     myEvent.nMCtrack=k+1;
                  }
               }
            } // end loop over stable particles
         }
      }
      // define initial state particle four-vectors
      double ee=*eBeamE;
      double pe= sqrt((ee+ME)*(ee-ME));
#ifdef CORRECT_FOR_TILT
      double pxe= - *beamtiltx0 *pe;
      double pye= - *beamtilty0 *pe;
#else
      double pxe= 0.;
      double pye= 0.;
#endif
      double pze = - sqrt(pe*pe-pxe*pxe-pye*pye);

      double ep=*eBeamP;
      static double const MP=0.9382720813;
      double pp= sqrt((ep+MP)*(ep-MP));
#ifdef CORRECT_FOR_TILT
      double pxp= *beamtiltx0 *pp;
      double pyp= *beamtilty0 *pp;
#else
      double pxp= 0.;
      double pyp= 0.;
#endif
      double pzp= sqrt(pp*pp-pxp*pxp-pyp*pyp);

      myEvent.eProtonBeam=*eBeamP;
      myEvent.eElectronBeam=*eBeamE;

      TLorentzVector ebeam_REC_lab(pxe,pye,pze,ee);
      TLorentzVector pbeam_REC_lab(pxp,pyp,pzp,ep);

      if(print) {
         cout<<"HERA beam energies: "<<ee<<" "<<ep<<" beam tilt: "<<*beamtiltx0<<" "<<*beamtilty0<<"\n";
         /* cout<<"Beam proton beam 4-vector\n";
         pbeam_REC_lab.Print();
         cout<<"Beam electron beam 4-vector\n";
         ebeam_REC_lab.Print(); */
      }

      // check HV conditions
      // not yet

      // trigger information
      double prob_rw=1.0,prob_ac=1.0;
      H1TrigInfo *trigInfo=dynamic_cast<H1TrigInfo *>
         (H1DBManager::Instance()->GetDBEntry(H1TrigInfo::Class()));
      // save all prescales etc for this run
      if(!trigInfo) {
         cout<<"TrigInfo not found!!!\n";
      }
      const Int_t *prescales=trigInfo->GetPrescales();
      const Int_t *enabled=trigInfo->GetEnabledSubTriggers();
      for(int i=0;i<4;i++) {
         myEvent.l1l2l3ac[i]=0;
         myEvent.l1l2l3rw[i]=0;
         for(int j=0;j<32;j++) {
            int st=i*32+j;
            if(!enabled[st]) continue;
            if(l1l2l3ac[st]) {
               myEvent.l1l2l3ac[i]|=(1<<j);
               prob_ac *= (1.-1./prescales[st]);
            }
            if(l1l2rw[st] && l1l3rw[st]) {
               myEvent.l1l2l3rw[i]|=(1<<j);
               prob_rw *= (1.-1./prescales[st]);
            }
         }
      }
      myEvent.trigWeightRW=(prob_rw<1.0) ? (1./(1.-prob_rw)) : 0.0;
      myEvent.trigWeightAC=(prob_ac<1.0) ? (1./(1.-prob_ac)) : 0.0;

      // background and noise finders
      myEvent.ibg=*ibg;
      myEvent.ibgfm=*ibgfm;
      myEvent.ibgam=*ibgam;
      myEvent.iqn=*iqn;

      // find primary vertex
      TVector3 beamSpot(*beamx0,*beamy0,0.);
      bool havePrimaryVertex=false;
      TVector3 primaryVertex(beamSpot);
      bool haveSimulatedVertex=false;
      TVector3 simulatedVertex(beamSpot);
      for(int i=0;i<vertex.GetEntries();i++) {
         Int_t type=vertex[i]->GetVertexType();
         if(type==0) {
            havePrimaryVertex=true;
            primaryVertex=vertex[i]->GetPosition();
         } else if(type==1) {
            haveSimulatedVertex=true;
            simulatedVertex=vertex[i]->GetPosition();
         }
      }
      myEvent.vertexType = *ivtyp;
      if(havePrimaryVertex) {
         myEvent.vertex[0]=primaryVertex.X();
         myEvent.vertex[1]=primaryVertex.Y();
         myEvent.vertex[2]=primaryVertex.Z();
      } else {
         myEvent.vertex[0]=-999.;
         myEvent.vertex[1]=-999.;
         myEvent.vertex[2]=-999.;
      }
      if(haveSimulatedVertex) {
         myEvent.simvertex[0]=simulatedVertex.X();
         myEvent.simvertex[1]=simulatedVertex.Y();
         myEvent.simvertex[2]=simulatedVertex.Z();
      } else {
         myEvent.simvertex[0]=-999.;
         myEvent.simvertex[1]=-999.;
         myEvent.simvertex[2]=-999.;
      }
      myEvent.beamSpot[0]=beamSpot.X();
      myEvent.beamSpot[1]=beamSpot.Y();
      myEvent.beamTilt[0]=*beamtiltx0;
      myEvent.beamTilt[1]=*beamtilty0;
      if(print) {
         cout<<"ivtyp="<<*ivtyp
             <<" beam spot: "<<beamSpot.X()<<" "<<beamSpot.Y();
         if(havePrimaryVertex) {
            cout<<" reconstructed primary vertex:"
                <<" "<<primaryVertex.X()
                <<" "<<primaryVertex.Y()
                <<" "<<primaryVertex.Z();
         }
         if(haveSimulatedVertex) {
             cout<<" simulated vertex:"
                 <<" "<<simulatedVertex.X()
                 <<" "<<simulatedVertex.Y()
                 <<" "<<simulatedVertex.Z();
         }
         cout<<"\n";
         cout<<"number of part cand: "<<partCand.GetEntries()<<"\n";
      }

      // find scattered electron as identified EM particle with highest PT in SpaCal
      bool haveScatteredElectron=false;
      TLorentzVector escat0_REC_lab;
      int scatteredElectron=-1;
      double ptMax=0;
      for(int i=0;i<partCand.GetEntries();i++) {
         H1PartCand *cand=partCand[i];
         H1PartEm const *elec=cand->GetIDElec();
         if(elec && (elec->GetType()==4)) {
            TLorentzVector p= elec->GetFourVector();
            if(p.Pt()>ptMax) {
               escat0_REC_lab = p;
               scatteredElectron=i;
               haveScatteredElectron=true;
            }
         }
      }
      
      // add EM particles and neutrals in a cone around the electron
      TLorentzVector escatPhot_REC_lab(escat0_REC_lab);
      set<int> isElectron;
      if(scatteredElectron>=0) {
         isElectron.insert(scatteredElectron);
         for(int i=0;i<partCand.GetEntries();i++) {
            if(i==scatteredElectron) continue;
            H1PartCand *cand=partCand[i];
            H1PartEm const *elec=cand->GetIDElec();
            if(elec) {
               TLorentzVector p= elec->GetFourVector();
               if(p.DeltaR(escat0_REC_lab)<ELEC_ISOLATION_CONE) {
                  escatPhot_REC_lab += p;
                  isElectron.insert(i);
               }
            } else if(!cand->GetTrack()) {
               TLorentzVector p= cand->GetFourVector();
               if(p.DeltaR(escat0_REC_lab)<ELEC_ISOLATION_CONE) {
                  escatPhot_REC_lab += p;
                  isElectron.insert(i);
               }
            }
         }
      }

      myEvent.elecEradREC=escatPhot_REC_lab.E()-escat0_REC_lab.E();
      myEvent.elecPxREC=escatPhot_REC_lab.X();
      myEvent.elecPyREC=escatPhot_REC_lab.Y();
      myEvent.elecPzREC=escatPhot_REC_lab.Z();
      myEvent.elecEREC=escatPhot_REC_lab.E();

      GetKinematics(ebeam_REC_lab,pbeam_REC_lab,escatPhot_REC_lab,
                    &myEvent.xREC,&myEvent.yREC,&myEvent.Q2REC);

      TLorentzRotation boost_REC_HCM=BoostToHCM(ebeam_REC_lab,pbeam_REC_lab,escatPhot_REC_lab);
      TLorentzVector q_REC_lab(ebeam_REC_lab-escatPhot_REC_lab);

      // calculate inclusive HFS 4-vector and track selection
      // exclude particles counted as electron
      TLorentzVector hfs;
      myEvent.nRECtrackAll=0;
      myEvent.nRECtrack=0;

      myEvent.nRECfstFitted=fstFittedTrack.GetEntries();
      myEvent.nRECfstSelected=0;

      vector<int> trackType(10);
      int nPart=partCand.GetEntries();

      nPart += fstFittedTrack.GetEntries();

      for(int i=0;i<nPart;i++) {
         H1PartCand *cand=0;
         H1FSTFittedTrack *fstTrack=0;
         TLorentzVector p;
         if(i<partCand.GetEntries()) {
            cand=partCand[i];
            p=cand->GetFourVector();
         } else {
            fstTrack=fstFittedTrack[i-partCand.GetEntries()];
            p=fstTrack->GetFourVector(M_CHARGED_PION);
         }
         // ignore particles counted with scattered electron
         if(cand && isElectron.find(i)!=isElectron.end()) continue;

         // exclude particles close to electron
         if(haveScatteredElectron &&
            (p.DeltaR(escat0_REC_lab)<ELEC_ISOLATION_CONE)) continue;

         if(cand) {
            // only particle candidates belong to the calibrated HFS
            hfs += p;
         }

         H1PartSelTrack const *track=0;
         if(cand) track=cand->GetIDTrack();
         if(track || fstTrack) {
            if(haveScatteredElectron) {
               TLorentzVector h=p;
               if(track) h=track->GetFourVector();
               double log10z=TMath::Log10((h*pbeam_REC_lab)/(q_REC_lab*pbeam_REC_lab));
               // boost to hadronic-centre-of-mass frame
               TLorentzVector hStar = boost_REC_HCM*h;
               double etaStar=hStar.Eta();
               double ptStar=hStar.Pt();
               double phiStar=hStar.Phi();
               int type=0;
               int charge=0;
               if(track) {
                  if(track->IsCentralTrk()) type =1;
                  else if(track->IsCombinedTrk()) type=2;
                  else if(track->IsForwardTrk()) type =3;
                  else if(track->IsFSTTrk()) type=4;
                  else if(track->IsBSTTrk()) type =5;
                  charge=track->GetCharge();
               } else if(fstTrack) {
                  // do some track selection here
                  // (1) tracks shall be a primary track
                  H1Vertex const *v=fstTrack->GetVertex();
                  if(floatEqual(v->X(),myEvent.vertex[0])&&
                     floatEqual(v->Y(),myEvent.vertex[1])&&
                     floatEqual(v->Z(),myEvent.vertex[2])) {
                     type=4;
                  }
                  // (2) minimum transverse momentum of 0.1 GeV
                  if(fstTrack->GetPt()<0.1) {
                     type=0;
                  }
                  // (3) momentum vector shall be incompatible with 
                  //  any other central, combined or forward track
                  if(type) {
                     charge=fstTrack->GetCharge();
                     TVector3 p1=fstTrack->GetMomentum();
                     TMatrix V1=fstTrack->GetMomentumCovar();
                     for(int j=0;j<partCand.GetEntries();j++) {
                         H1PartCand *candJ=partCand[j];
                         H1PartSelTrack const *selTrackJ=candJ->GetIDTrack();
                         H1PartCand const *partCandJ=
                            selTrackJ ? (selTrackJ->GetParticle()) : 0;
                         H1Track const *trackJ=partCandJ ? partCandJ->GetTrack() : 0;
                         if(trackJ) {
                            TVector3 p2=trackJ->GetMomentum();
                            TMatrix V2=trackJ->GetMomentumCovar();
                            TMatrixD sum(V1+V2);
                            TMatrixD Vinv(TMatrixD::kInverted,V1+V2);
                            TVector3 d(p1-p2);
                            double chi2=d.Dot(Vinv*d);
                            //if(print) cout<<i<<" "<<j<<" "<<chi2;
                            if(chi2<30.) {
                               //if(print) cout<<" [reject]";
                               type=0;
                            }
                            //if(print) cout<<"\n";
                         }
                     }
                  }
                  if(type) {
                     myEvent.nRECfstSelected++;
                  }
               }
               trackType[type]++;
               if(type && (myEvent.nRECtrack<MyEvent::nRECtrack_MAX)) {
                  if(print) {
                     cout<<i<<" Track "<<myEvent.nRECtrackAll
                         <<" "<<charge*type
                         <<" etaLab="<<h.Eta()
                         <<" ptLab="<<h.Pt()
                         <<" phiLab="<<h.Phi()
                         <<" ptStar="<<ptStar
                         <<" etaStar="<<etaStar
                         <<" phiStar="<<phiStar
                         <<" log10(z)="<<log10z
                         <<"\n";
                  }
                  myEvent.nRECtrackAll++;
                  int k=myEvent.nRECtrack;
                  myEvent.typeChgREC[k]=charge*type;
#ifdef SAVE_LAB_MOMENTA
                  myEvent.pxREC[k]=h.X();
                  myEvent.pyREC[k]=h.Y();
                  myEvent.pzREC[k]=h.Z();
#endif
                  myEvent.ptStarREC[k]=hStar.Pt();
                  myEvent.etaStarREC[k]=hStar.Eta();
                  myEvent.phiStarREC[k]=hStar.Phi();
                  myEvent.log10zREC[k]=log10z;
                  myEvent.nucliaREC[k]=1.;
                  myEvent.momREC[k]=h.Vect();
                  myEvent.covREC[k].ResizeTo(3,3);
                  myEvent.imatchREC[k]=-2;
                  myEvent.dmatchREC[k]=-1.;
                  if(fstTrack) {
                     myEvent.covREC[k]=fstTrack->GetMomentumCovar();
                     myEvent.imatchREC[k]=-1;
                  } else {
                     H1PartCand const *partCandI=track->GetParticle();
                     H1Track const *trackI=partCandI ? partCandI->GetTrack():0;
                     if(trackI) {
                        myEvent.covREC[k]=trackI->GetMomentumCovar();
                        myEvent.imatchREC[k]=-1;
                     }
                  }
                  myEvent.nRECtrack=k+1;
               }
            }
         }
      }
      // match MC particles and REC particles
      // (1) for each REC particle, find the best MC particle
      //    [may result in multiple REC particles matched to the same MCpart]
      //    matching is perfomed by selecting the MC particle which 
      //    gives the lowest chi**2 when comparing the momenta
      //
      //  this sets:  myEvent.dmatchREC[]  -> lowest chi**2
      //              myEvent.imatchREC[]  -> best matching MC particle
      for(int iREC=0;iREC<myEvent.nRECtrack;iREC++) {
         // skip tracke where momentum covariance is not known
         //  -> these will never be matched
         if(myEvent.imatchREC[iREC]!=-1) continue;
         TMatrixDSym Vsym(3);
         for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
               Vsym(i,j)=myEvent.covREC[iREC](i,j);
            }
         }
         TMatrixDSymEigen ODO(Vsym);
         TVectorD ev=ODO.GetEigenValues();
         //if(print) ev.Print();
         TMatrixD O=ODO.GetEigenVectors();
         TMatrixD Ot(TMatrixD::kTransposed,O);

         //TMatrixD Vinv(TMatrixD::kInverted,myEvent.covREC[iREC]);
         for(int jMC=0;jMC< myEvent.nMCtrack;jMC++) {
            TVector3 d=myEvent.momREC[iREC]- myEvent.partMC[jMC]->GetMomentum();
            TVector3 Otd=Ot*d;
            double chi2=0.;
            for(int i=0;i<3;i++) {
               if(ev[i]>=1.E-6*ev[0]) {
                  chi2+=Otd[i]*Otd[i]/ev[i];
               }
            }
            //double chi2simple=d.Dot(Vinv*d);
            //if(print) cout<<iREC<<" "<<jMC<<" "<<chi2<<" "<<chi2simple<<"\n";
            if((jMC==0)||(chi2<myEvent.dmatchREC[iREC])) {
               myEvent.dmatchREC[iREC]=chi2;
               myEvent.imatchREC[iREC]=jMC;
            }
         }
      }
      // (2) set pointers from MC to REC
      //     and remove duplicates
      for(int iREC=0;iREC<myEvent.nRECtrack;iREC++) {
         int iMC=myEvent.imatchREC[iREC];
         if(iREC<0) continue; // no match for this particle
         int jREC=myEvent.imatchMC[iMC];
         if(jREC>=0) {
            // duplicate match for this particle
            // iREC and jREC both are pointing to the same particle
            // compare matching distance
            if(myEvent.dmatchREC[jREC]<myEvent.dmatchREC[iREC]) {
               // old match is better
               // invalidate pointer REC->MC
               myEvent.imatchREC[iREC]=-2-iMC;
            } else {
               // new match is better
               // invalidate old pointer REC->MC
               myEvent.imatchREC[jREC]=-2-iMC;
               // save new pointer MC->REC
               myEvent.imatchMC[iMC]=iREC;
            }
         } else {
            // save this match
            myEvent.imatchMC[iMC]=iREC;
         }
      }
      // now:
      //    myEvent.imatchMC[] points to the best matching REC particle
      //                         <0 -> inefficiency
      //    myEvent.imatchREC[] points to the best matching MC particle
      //                         <0 -> fake track

      // for matched particles, calculate extra weight for
      // nuclear interaction probability correction
      for(int iREC=0;iREC<myEvent.nRECtrack;iREC++) {
         int iMC=myEvent.imatchREC[iREC];
         int part=0;
         if(iMC>=0) {
            int pdg=myEvent.idMC[iMC];
            if(pdg<0) pdg= -pdg;
            part= (pdg==211) ? 1 : ((pdg==321) ? 2 : 0);
         }
         if(part) {
            myEvent.nucliaREC[iREC]=
               H1NuclIACor::GetWeight
               (part,myEvent.typeChgREC[iREC]>0 ? 1 : -1,
                myEvent.momREC[iREC].Pt(),
                myEvent.momREC[iREC].Phi(),myEvent.momREC[iREC].Theta(),
                0.0 /* dca */,H1SelVertex::GetPrimaryVertex()->Z());
         }
      }

      if(print) {
         for(int iREC=0;iREC<myEvent.nRECtrack;iREC++) {
            if(myEvent.imatchREC[iREC]>=0) {
               cout<<"REC track "<<iREC<<" is matched to MC particle "
                   << myEvent.imatchREC[iREC]<<" dmatch="
                   <<myEvent.dmatchREC[iREC]
                   <<" nuclIA="<<myEvent.nucliaREC[iREC]
                   <<"\n";
            }
         }
         for(int iREC=0;iREC<myEvent.nRECtrack;iREC++) {
            if(myEvent.imatchREC[iREC]<0) {
               cout<<"REC track "<<iREC<<" NOT matched "
                   << myEvent.imatchREC[iREC]<<" dmatch="
                   <<myEvent.dmatchREC[iREC]
                   <<" nuclIA="<<myEvent.nucliaREC[iREC]
                   <<"\n";
            }
         }
      }


      myEvent.hfsPxREC=hfs.X();
      myEvent.hfsPyREC=hfs.Y();
      myEvent.hfsPzREC=hfs.Z();
      myEvent.hfsEREC=hfs.E();

      if(haveScatteredElectron && print) {
         cout<<"reconstructed electron w/o photons in lab: ";
         escat0_REC_lab.Print();
         cout<<"reconstructed electron with photons in lab: ";
         escatPhot_REC_lab.Print();
      }

      if(print) {
         for(size_t type=0;type<trackType.size();type++) {
            cout<<" "<<trackType[type];
         }
         cout<<"\n";
      }

      if(print) {
         print--;
      }

      output->Fill();
   }



    // Summary
    cout << "\nProcessed " << eventCounter
        << " events\n\n";


    // Write histogram to file
    output->Write();
    file.Close();

    cout << "Histograms written to " << opts.GetOutput() << endl;

    return 0;
}
