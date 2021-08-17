#include "FidVolCut.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <TMath.h>

#include "H1Steering/H1SteerManager.h"
#include "H1Calculator/H1Calculator.h"
#include "H1Calculator/H1CalcElec.h"
#include "H1Calculator/H1CalcEvent.h"

using std::cout;
using std::endl;
using std::setw;

ClassImp(FidVolCut)

//________________________________________________________
//
// Class FidVolCut - cuts out the fiducial volume
// where the trigger efficiency is less than 100 %
//________________________________________________________
//

FidVolCut::FidVolCut() : H1Cut()
{
}

//-------------------------------------------------------
FidVolCut::FidVolCut(TString name) : H1Cut()
{
    Int_t NumPeriods;
    Int_t*   tempi;
    Float_t* tempf;

    if(gH1SteerManager != 0) {
        gH1SteerManager->SetFilename("FidVolCut.steer");

    	// choose the fiducial volume cut depending on the run period
    	if (name.Contains("HERA1")) { // HERA I 
    	  cout << "Reading fiducial volume cut for period 99/00. Name = " << name.Data() << endl;
    	  gH1SteerManager->SetClassname("FidVolCut", "FidVolCut9900");
    	}
        else if (name.Contains("HERA2")) {  // HERA II
    	  cout << "Reading fiducial volume cut for period 03/07. Name = " << name.Data() << endl;
    	  gH1SteerManager->SetClassname("FidVolCut", "FidVolCut0307");
    	}
        else {
    	   Error("FidVolCut", "Can not initialize fiducial volume cuts for unknown period. Name is %s", name.Data());
    	}

        NumPeriods = gH1SteerManager->ReadInt("fNumPeriods");

        tempi = new Int_t[NumPeriods];
        tempf = new Float_t[NumPeriods];
      
        gH1SteerManager->ReadArray("fRunStart", tempi);
        for(int i = 0 ; i < NumPeriods; ++i) fRunStart += tempi[i];
        gH1SteerManager->ReadArray("fRunEnd", tempi);
        for(int i = 0 ; i < NumPeriods; ++i) fRunEnd += tempi[i];
        gH1SteerManager->ReadArray("fPhiStart", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fPhiStart += tempf[i];
        gH1SteerManager->ReadArray("fPhiEnd", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fPhiEnd += tempf[i];
        gH1SteerManager->ReadArray("fZImpStart", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fZImpStart += tempf[i];
        gH1SteerManager->ReadArray("fZImpEnd", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fZImpEnd += tempf[i];
        gH1SteerManager->ReadArray("fElecEStart", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fElecEStart += tempf[i];
        gH1SteerManager->ReadArray("fElecEEnd", tempf);
        for(int i = 0 ; i < NumPeriods; ++i) fElecEEnd += tempf[i];
    }
    else {
        cout << "no steer file determined" << endl;
    }
}
//-------------------------------------------------------

//-------------------------------------------------------
FidVolCut::FidVolCut(const FidVolCut& BluePrint) : H1Cut(BluePrint)
{
    fRunStart   = BluePrint.fRunStart;
    fRunEnd     = BluePrint.fRunEnd;
    fPhiStart   = BluePrint.fPhiStart;
    fPhiEnd     = BluePrint.fPhiEnd;
    fZImpStart  = BluePrint.fZImpStart;
    fZImpEnd    = BluePrint.fZImpEnd;
    fElecEStart = BluePrint.fElecEStart;
    fElecEEnd   = BluePrint.fElecEEnd;
}
//-------------------------------------------------------

//-------------------------------------------------------
FidVolCut::~FidVolCut()
{
}
//-------------------------------------------------------


//-------------------------------------------------------
Bool_t FidVolCut::PassesCut(Int_t index) const
{
    // The fiducial volume cut

    return FiducialVolumeCut();
}
//-------------------------------------------------------


//-------------------------------------------------------  
Bool_t FidVolCut::FiducialVolumeCut() const
{
    // Removes areas of the EM Calo where triggering efficiency is not at 100%
    // Returns fTRUE if all the FV criterea listed in fiducial.txt file are met; if 
    // there is no electron; or if running over a Monte Carlo.
    // fiducial .txt is written: RunMin RunMax PhiMin PhiMax ZimpMin ZimpMax EMin EMax

    Bool_t flag = kTRUE;

    Float_t Zimpact = gH1Calc->Elec()->GetZimpact();
    //  Float_t Phi = gH1Calc->Elec()->GetFirstElectron().Phi();
    Float_t Phi = gH1Calc->Elec()->GetPhiAtLarEdge();
    Float_t E = gH1Calc->Elec()->GetFirstElectron().E();
    Float_t RunNumber = 0;
  
    if(!(gH1Calc->IsMC())){
        RunNumber  = gH1Calc->GetRunNumber();
    } else {
        RunNumber = gH1Calc->Event()->GetMCRunNumber();// leave out for now
    }

    for(int i = 0; i < fRunStart.GetEntries(); ++i){
        if(RunNumber < fRunStart[i]) continue;
        if(RunNumber > fRunEnd[i]) continue;
        if(Phi*( 180/TMath::Pi() ) < fPhiStart[i]) continue;
        if(Phi*( 180/TMath::Pi() ) > fPhiEnd[i]) continue;
        if(Zimpact < fZImpStart[i]) continue;
        if(Zimpact > fZImpEnd[i]) continue;
        if(E < fElecEStart[i]) continue;
        if(E > fElecEEnd[i]) continue;
        flag = kFALSE;
    }

    return flag ? true : BooleanCut(flag);
}
//-------------------------------------------------------


//-------------------------------------------------------
void FidVolCut::Print(ostream* os) const
{
    // Print out details of cuts implemented
    // Name, and any ranges involved

    *os << "The Fiducial Volume Cut used the ranges :" << endl;
    *os << setw(11) << "Run Start"   << setw(11) << "Run End"
        << setw(11) << "Phi Start"   << setw(11) << "Phi End" 
        << setw(11) << "Zimp Start"  << setw(11) << "Zimp End"
        << setw(12) << "ElecE Start" << setw(12) << "ElecE End"
        << endl;

    for(int i = 0 ; i < fRunStart.GetEntries(); ++i){
        *os << setw(11) << fRunStart[i]   << setw(11) << fRunEnd[i] 
            << setw(11) << fPhiStart[i]   << setw(11) << fPhiEnd[i]
            << setw(11) << fZImpStart[i]  << setw(11) << fZImpEnd[i]
            << setw(12) << fElecEStart[i] << setw(12) << fElecEEnd[i]
            << endl;
    }
}
//-------------------------------------------------------

void FidVolCut::PrintDebug() const 
{
  // print debug info - minimalist for now
  cout << "FidVolCut: Failed "<<endl;
}