#include <TLorentzVector.h>

class H1PartMCArrayPtr;

class TDetectQedc {
 public:
   // create this class to detect QEDC events
   TDetectQedc(H1PartMCArrayPtr &mcpart);
   void Print(H1PartMCArrayPtr &mcpart);
   Bool_t IsElectronFound(void) const;
   Bool_t IsPhotonQedcFound(void) const;
   Int_t GetNumISR(void) const { return fNisr; }
   Int_t GetNumFSR(void) const { return fNfsr; }
   Bool_t IsQedcEvent(void) const;
   Int_t GetCutFlags(Int_t i) const { return fCutsQedc[i]; }
   TLorentzVector const &GetElectronNoFSR(void) const { return fElectronNoFSR; }
   TLorentzVector const &GetElectron(void) const { return fElectron; }
   TLorentzVector const &GetPhoton(int i) const { return fPhoton[i]; }
   TLorentzVector const &GetW(void) const { return fW; }
   TLorentzVector const &GetInvis(void) const { return fInvis; }

   Double_t GetMstring(void) const { return fMstring; }
   Double_t GetMcluster(void) const { return fMcluster; }
   Int_t GetNumString(void) const { return fNstring; }
   Int_t GetNumCluster(void) const { return fNcluster; }
   Int_t GetMesonPdg(void) const { return fMesonPdg; }
   TLorentzVector const &GetMeson(void) const { return fMeson; }
   // photon types
   enum PHOTON_TYPE { PHOTON_QEDC=0, PHOTON_ISR=1, PHOTON_FSR=2,
                      PHOTON_NTYPE=3};
 protected:
   // check cuts for a given electron, photon pair
   Int_t GetQedcCuts(TLorentzVector const &electron,
                    TLorentzVector const &photon) const;
   // QEDC generator cuts
   Float_t fELIE,fAGMA,fAGMI,fELIG,fPTMA,fVMIN,fWMIN,fWMAX,fACO;
   static Float_t const kELIE_DEFAULT,kAGMA_DEFAULT,kAGMI_DEFAULT,
      kELIG_DEFAULT,kPTMA_DEFAULT,kVMIN_DEFAULT,kWMIN_DEFAULT,
      kWMAX_DEFAULT,kACO_DEFAULT;

   // index of electron and photon originating from QEDC
   // (only if identified as such)
   Int_t fIelectron,fIphotonQedc;
   // index of electron prior to FSR (COMPTON23)
   Int_t fIelectronNoFSR;
   // diffractive variables, to detect poor-mans pomeron
   Int_t fNcluster,fNstring;
   // detect diffractive (vector)-meson production
   Int_t fMesonPdg;
   TLorentzVector fMeson;
   // scalar mass of all strings
   Double_t fMstring,fMcluster;
   Double_t fMX,fMY;
   // number of ISR/FSR photons
   Int_t fNisr,fNfsr;
   // four-vector of scattered electron
   TLorentzVector fElectron;
   TLorentzVector fElectronNoFSR;
   // four-vectors of photons
   TLorentzVector fPhoton[PHOTON_NTYPE];
   // cut flags for each electron-photon pairing
   // if cut is zero, then the e,gamma pair survived all QEDC generator cuts 
   Int_t fCutsQedc[PHOTON_NTYPE];
   // incoming particle 4-vector
   TLorentzVector fEbeam,fEPbeam,fW,fInvis;
};
