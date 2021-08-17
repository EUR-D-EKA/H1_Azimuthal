#ifndef __FIDVOLCUT__H
#define __FIDVOLCUT__H

#include "H1Cuts/H1Cut.h"
#include "H1Arrays/H1ArrayI.h"
#include "H1Arrays/H1ArrayD.h"

class FidVolCut : public H1Cut
{

    public:

        FidVolCut();
        explicit FidVolCut(TString name);
        FidVolCut(const FidVolCut& BluePrint);
        ~FidVolCut();

        Bool_t PassesCut(Int_t index = -1) const;

        void Print(ostream* os) const;  

        FidVolCut* Clone() const {return new FidVolCut(*this);}

        Bool_t FiducialVolumeCut() const;

        void SetIndex(const Int_t Index = -1){}

        void PrintDebug() const;

    private:

        H1ArrayI fRunStart;     //!
        H1ArrayI fRunEnd;       //!
        H1ArrayD fPhiStart;     //!
        H1ArrayD fPhiEnd;       //!
        H1ArrayD fZImpStart;    //!
        H1ArrayD fZImpEnd;      //!
        H1ArrayD fElecEStart;   //!
        H1ArrayD fElecEEnd;     //!
  
    private:
  
        ClassDef(FidVolCut, 1) // Cut class containing cuts on more 
                               // complicated variables (non-generic access)
};
#endif