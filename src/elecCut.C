/****************************************************************
 * File: elecCut.C 
 *
 * Provides easy access to some more hopefully sophisticated cuts 
 * for the detection of the scattered electron. 
 ****************************************************************/ 

#include <iostream> 
#include <sstream>
#include <stdlib.h> 

//root Header
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h> 
#include <TF1.h>
#include <TGraph.h> 
#include <TGraphErrors.h> 
#include "TGraphAsymmErrors.h"
#include <TLatex.h>
#include <TMath.h> 
#include <TPostScript.h>
#include <TLine.h> 

#include "TAttLine.h" 
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TAttText.h"
#include "TClass.h"

#include "TMinuit.h" 
#include "TMatrixD.h"
#include "TVectorF.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TVirtualFitter.h"
#include "TEllipse.h"
#include "TPave.h"

#include "H1Clusters/H1Cell.h"
#include "H1Geom/H1CellGeometry.h"
#include "H1Geom/H1SpaCalCellGeo.h"
#include "H1Geom/H1CaloGeometry.h"
#include "H1Geom/H1SpacalGeometry.h"
#include "elecCut.h"

static H1FloatPtr ElecTheta("ElecTheta");      // theta of scattered electron

//ClassImp(elecCut) 

/***************************************************************
 * initElecCut --  
 * for MC a run number must be provided in order to set the 
 * run range of the mc. 
 */
elecCut::elecCut(int RunNumber_) {
  dbg=kTRUE;  
  nDeadCells = 0; 
  nDeadTrigCells = 0;
  nCircularCuts=0; 

  RunNumber=0; 
  RunRange=0; 

 //std circle (always present): 
  if(dbg) printf("elecCut::elecCut Initialising std. circular cut\n"); 
  nCircularCuts = 1; 
  circCutsRadius[0] = 12 ;
  circCutsXzero[0] = -2.025; 
  circCutsYzero[0] = 0; 
  
  if(RunNumber_ !=0) { 
    printf("elecCut::elecCut: Initialising Cuts for given runNumber=%d\n", RunNumber); 
    initElecCut(RunNumber_);
  }
  dbg=kFALSE; 
} 

/*****************************************************************
 * initElecCut -- initialise the box cuts according to the given 
 * run number. At the moment three different run periods are used
 * which have different initialsing functions. These function are 
 * from here according to the given run number. 
 ****************************************************************/
int elecCut::initElecCut(int RunNumber_) { 
  int newRunRange=0; 

   if(RunNumber_!=RunNumber) {
     if(dbg) printf(" elecCut::initElecCut: RunNumber has changed from %d to %d check if SpaCal-Cuts needs to be changed\n", RunNumber, RunNumber_); 
    
     if(RunNumber_ >=367257 && RunNumber_ <= 392213)  //2004     (e+ running)  (2004 until end: 398818)
       newRunRange=2004 ;
     else if(RunNumber_ >398286 && RunNumber_ <= 436893) //2004-2005 (e- running)
       newRunRange=2005;
     else if(RunNumber_ >436893 && RunNumber_ <= 492558 )  //2006: 466977
       newRunRange=2006;
     else if(RunNumber_ >492558 && RunNumber_ <=500611  )  {//2007:
     	newRunRange=2007;
//      	printf("RunRange 2006b-2007 requested. Not analysed yet. Using Spacal Cuts of 2006 instead.\n"); 
     }
     else { 
       printf(" elecCut::initElecCut: -- WARNING -- Run %d is in an unknown RunRange. Using current RunRange %d instead\n",RunNumber_, RunRange); 
       newRunRange = RunRange;
     }

     if(newRunRange!=RunRange) { 
       RunNumber = RunNumber_; 
       RunRange = newRunRange; 
       if(dbg) printf(" elecCut::initElecCut: RunNumber has changed and RunRange also -> Initialising new Cuts for the SpaCal\n"); 
       if(dbg) printf(" elecCut::initElecCut: Initilising SpaCalCut for the year %d\n", newRunRange);
       clearDeadCellPositions(); 
       if(newRunRange == 2004)  initCuts2004(); 
       else if(newRunRange == 2005)  initCuts2005(); 
       else if(newRunRange == 2006)  initCuts2006();
       else if(newRunRange == 2007)  initCuts2007();
       else printf("elecCut::initElecCut: -- WARNING -- RunRange %d is not known. NO CUTS initialised\n", newRunRange);
       RunRange = newRunRange;
     }
     RunNumber = RunNumber_; 
   }
   //}
  return 1; 
} 


/*****************************************************************
 * initCuts2004 -- box cuts for the year 2004
 *
 *
 ****************************************************************/
void elecCut::initCuts2004() {

  printf("elecCut::elecCut: Initialising deadSpacal and TriggerCells for year 2004\n");
 
  clearDeadCellPositions();

  nDeadCells = 6; 
  nDeadTrigCells = 13; 
  
  deadCells[0] = 487;
  deadCells[1] = 128; 
  deadCells[2] = 403;
  deadCells[3] = 197;
  deadCells[4] = 306;
  deadCells[5] = 763;

  //initialise array of dead trigger cells: 
  deadTrigCells[0] = 22; 
  deadTrigCells[1] = 23; 
  deadTrigCells[2] = 45;
  deadTrigCells[3] = 46; 
  deadTrigCells[4] = 59; 
  deadTrigCells[5] = 60; 
  deadTrigCells[6] = 66; 
  deadTrigCells[7] = 67; 
  deadTrigCells[8] = 94; 
  deadTrigCells[9] = 95; 
  deadTrigCells[10] = 103; 
  deadTrigCells[11] = 104; 
  deadTrigCells[12] = 89;  

  //circular cuts: 
  printf("Implementing an additional circ. cut (before ncuts=%d) \n", nCircularCuts); 
  nCircularCuts++; 
  circCutsRadius[nCircularCuts-1] = 10.4;
  circCutsXzero[nCircularCuts-1] = 0.260701;
  circCutsYzero[nCircularCuts-1] = 3.15767; 

  nCircularCuts++; 
  circCutsRadius[nCircularCuts-1] = 7.87295; 
  circCutsXzero[nCircularCuts-1] = 1.63656;
  circCutsYzero[nCircularCuts-1] = -4.20146; 

  printf("There are now %d circ. cuts implemented\n", nCircularCuts); 

  H1DBManager::Instance()->StartRun(367258);
  fillDeadCellPositions();
  RunRange = 2004; 
  return; 
} 

/*****************************************************************
 * initCuts2005 -- box cuts for the year 2005 
 *
 *
 ****************************************************************/
void elecCut::initCuts2005(){ 
  dbg=kFALSE; 

  clearDeadCellPositions(); 

   printf("elecCut::elecCut: Initialising deadSpacal and TriggerCells for year 2005\n");
  
  nDeadCells = 11; 
  nDeadTrigCells = 17; 
  
  //initialise array of dead cells: 
  deadCells[0] = 75; 
  deadCells[1] = 76; 
  deadCells[2] = 223; 
  deadCells[3] = 306; 
  deadCells[4] = 190; 
  deadCells[5] = 258; 
  deadCells[6] = 403;
  deadCells[7] = 371; 
  deadCells[8] = 1010;
  deadCells[9] = 915;
  deadCells[10] = 789;

  deadTrigCells[0] = 22; 
  deadTrigCells[1] = 23; 
  deadTrigCells[2] = 45; 
  deadTrigCells[3] = 46; 
  deadTrigCells[4] = 44; 
  deadTrigCells[5] = 59;
  deadTrigCells[6] = 60;
  deadTrigCells[7] = 94;
  deadTrigCells[8] = 95;

  deadTrigCells[9] = 590;
  deadTrigCells[10] =591;
  deadTrigCells[11] = 691;
  deadTrigCells[12] = 692;
  deadTrigCells[13] = 763;

  deadTrigCells[14] = 241;

  deadTrigCells[15] = 644;
  deadTrigCells[16] = 547;

   //circular cuts: 
  printf("Implementing an additional circ. cut (before ncuts=%d) \n", nCircularCuts); 
  nCircularCuts++; 
  circCutsRadius[nCircularCuts-1] = 10.2;
  circCutsXzero[nCircularCuts-1] = 0.258829;
  circCutsYzero[nCircularCuts-1] = 2.09418; 

  nCircularCuts++; 
  circCutsRadius[nCircularCuts-1] = 7.7; 
  circCutsXzero[nCircularCuts-1] = 2.46162;
  circCutsYzero[nCircularCuts-1] = -2.89197; 

  printf("There are now %d circ. cuts implemented\n", nCircularCuts); 


  H1DBManager::Instance()->StartRun(398820);
  fillDeadCellPositions();
  RunRange = 2005; 
  return; 
}

/*****************************************************************
 * initCuts2006 -- box cuts for the year 2006
 *
 *
 ****************************************************************/
void elecCut::initCuts2006() { 
 printf("elecCut::elecCut: Initialising deadSpacal and TriggerCells for year 2006\n");

  clearDeadCellPositions(); 

  nDeadCells = 2; 
  deadCells[0] = 403;
  deadCells[1] = 297;
  
  nDeadTrigCells = 4;
  deadTrigCells[0] = 22;
  deadTrigCells[1] = 23;
  deadTrigCells[2] = 45 ; 
  deadTrigCells[3] = 46 ; 

 
  //deadTrigCells[0] = 6; 
  //deadTrigCells[1] = 9; 
  //deadTrigCells[2] = 18; 
  //deadTrigCells[3] = 19;
  //deadTrigCells[6] = 26;
  //deadTrigCells[7] = 27;
  //deadTrigCells[10] = 21 ;
  //deadTrigCells[11] = 24 ;
  //deadTrigCells[10] = 59 ;
  //deadTrigCells[11] = 60 ;
  //deadTrigCells[14] = 76 ;
  //deadTrigCells[15] = 77 ;
  //deadTrigCells[16] = 94 ;
  //deadTrigCells[17] = 95 ;

   H1DBManager::Instance()->StartRun(436894);
  fillDeadCellPositions();


  //circular cuts:
  printf("Implementing an additional circ. cut (before ncuts=%d) \n", nCircularCuts);
  nCircularCuts++;
  circCutsRadius[nCircularCuts-1] = 10.4;
  circCutsXzero[nCircularCuts-1] = 0.337138;
  circCutsYzero[nCircularCuts-1] = 2.01292;

  nCircularCuts++;
  circCutsRadius[nCircularCuts-1] = 7.26146;
  circCutsXzero[nCircularCuts-1] = 2.47736;
  circCutsYzero[nCircularCuts-1] = -4.52853;

  printf("There are now %d circ. cuts implemented\n", nCircularCuts);
   RunRange = 2006;

  return;
}

/*****************************************************************
 * initCuts2007 -- box cuts for the year 2007
 * - upt to now 12.09.07 the code is a copy of
 * the one from 2006 !!!
 ****************************************************************/
void elecCut::initCuts2007() {
 printf("elecCut::elecCut: Initialising deadSpacal and TriggerCells for year 2007 - be aware this is just a copy of 2006-cuts! \n");

  clearDeadCellPositions(); 

  nDeadCells = 2; 
  deadCells[0] = 403;
  deadCells[1] = 297;
  
  nDeadTrigCells = 4;
  deadTrigCells[0] = 22;
  deadTrigCells[1] = 23;
  deadTrigCells[2] = 45 ; 
  deadTrigCells[3] = 46 ; 

 
  //deadTrigCells[0] = 6; 
  //deadTrigCells[1] = 9; 
  //deadTrigCells[2] = 18; 
  //deadTrigCells[3] = 19;
  //deadTrigCells[6] = 26;
  //deadTrigCells[7] = 27;
  //deadTrigCells[10] = 21 ;
  //deadTrigCells[11] = 24 ;
  //deadTrigCells[10] = 59 ;
  //deadTrigCells[11] = 60 ;
  //deadTrigCells[14] = 76 ;
  //deadTrigCells[15] = 77 ;
  //deadTrigCells[16] = 94 ;
  //deadTrigCells[17] = 95 ;

   H1DBManager::Instance()->StartRun(492558);
  fillDeadCellPositions();


  //circular cuts:
  printf("Implementing an additional circ. cut (before ncuts=%d) \n", nCircularCuts);
  nCircularCuts++;
  circCutsRadius[nCircularCuts-1] = 10.4;
  circCutsXzero[nCircularCuts-1] = 0.337138;
  circCutsYzero[nCircularCuts-1] = 2.01292;

  nCircularCuts++;
  circCutsRadius[nCircularCuts-1] = 7.26146;
  circCutsXzero[nCircularCuts-1] = 2.47736;
  circCutsYzero[nCircularCuts-1] = -4.52853;

  printf("There are now %d circ. cuts implemented\n", nCircularCuts);
   RunRange = 2007;

  return;
}

/*******************************************************************
 * clearDeadCellPositions -- clear arrays with cell positions
 *
 *
 *******************************************************************/
void elecCut::clearDeadCellPositions() {
  nDeadCells = 0;
  nDeadTrigCells = 0;

  for(int i=0;i<MAX_DEAD_CELLS;i++) { 
    deadCells[i] = 0; 
    posDeadCellsX[i] =0; 
    posDeadCellsY[i] =0; 
  } 
  for(int i=0;i<MAX_DEAD_TRIG_CELLS;i++) {
    deadTrigCells[i] = 0; 
    posDeadTrigCellsX[i] = 0;
    posDeadTrigCellsY[i] = 0; 
  }

 //std circle (always present): 
  nCircularCuts = 1; 
  circCutsRadius[0] = 12 ;
  circCutsXzero[0] = -2.025; 
  circCutsYzero[0] = 0;

  for(int i=1;i<MAX_CIRC_CUTS;i++) {
    circCutsRadius[i]=0;
    circCutsXzero[i]=0;
    circCutsYzero[i]=0;
  }

  return; 
} 

/******************************************************************
 * fillDeadCellPositions -- use the given cell number to fill the
 * arrays with the cut bounds in order to cut out these cells. Used 
 * for all run ranges. 
 *****************************************************************/
void elecCut::fillDeadCellPositions() { 
  if(dbg) printf("elecCut::fillDeadCellPositions: Getting Cell Positions from DB for Run %d\n", RunNumber); 

  //make instance of DB specifiing the run number: 
  H1DBManager::Instance()->StartRun(RunNumber); 
  
  if(dbg)printf("elecCut::fillDeadCellPositions: Trying to access DB for Geometry information\n"); 
  H1CellGeometry* fGeometry; 
  H1CaloGeometry * caloGeom = static_cast<H1CaloGeometry *> (H1DBManager::Instance()->GetDBEntry(H1CaloGeometry::Class()));
  if(!caloGeom){
    printf("elecCut::fillDeadCellPositions: GetDBGeometry: No H1CaloGeometry available.\n");
  }
   
  for(int i=0;i<nDeadCells;i++) { 
    
    if(dbg)printf("elecCut::fillDeadCellPositions: Getting SpaCal Position Information of caloCell: %d\n",deadCells[i] ); 
    fGeometry = caloGeom->GetCell( H1CaloGeometry::kSpaCal,deadCells[i]);
    if(fGeometry == NULL) { 
      printf("elecCut::elecCut: Warning request to exclude not existing SpaCal Cell! Might give strange results\n"); 
      continue; 
    } 
    posDeadCellsX[i] = fGeometry->GetX();
    posDeadCellsY[i] = fGeometry->GetY();
    if(dbg)printf("elecCut::fillDeadCellPositions: SpaCal-Cell: %d; x=%f y=%f\n", deadCells[i],  posDeadCellsX[i], posDeadCellsY[i]); 
  }
  
  for(int i=0;i<nDeadTrigCells;i++) { 

    if(dbg) printf("elecCut::fillDeadCellPositions: Dead TrigCells: Getting SpaCal Position Information of caloCell: %d\n",deadTrigCells[i]); 
    fGeometry = caloGeom->GetCell( H1CaloGeometry::kSpaCal,deadTrigCells[i]);


    if(fGeometry == NULL) { 
      printf("elecCut::fillDeadCellPositions: Dead TrigCells: Warning request to exclude not existing SpaCal Cell! Might give strange results. Array Nmb: %d, CellNmb: %d\n", i, deadTrigCells[i]); 
      continue; 
    } 
    posDeadTrigCellsX[i] = fGeometry->GetX(); 
    posDeadTrigCellsY[i] = fGeometry->GetY(); 
    if(dbg)printf("elecCut::fillDeadCellPositions: Dead TrigCells:SpaCal-Cell: %d; x=%f y=%f\n",deadTrigCells[i], posDeadTrigCellsX[i], posDeadTrigCellsY[i]); 
  }

  return; 
}

/****************************************************************
 * goodElec -- select a good electron asking for: 
 *
 * 1. not inside of spacial cuts 
 * 2. cut on hadronic energy fraction
 * 3. track poiting to cluster 
 */
int elecCut::goodElec( H1PartEm const *elec, Int_t RunNumber_) { 
 
  //check if cuts are still ok: 
  initElecCut(RunNumber_); 


 
  //elec whether h1oo finder found it: 
  if(!elec->IsScatElec() ) { 
    return 0; 
  } 
  
  //make spatial cuts: 
  if( elecCutSpatial( elec->GetXClus(),  elec->GetYClus() )!=1){ 
     return 0;
  } 
  
  //make cut on hadronic energy part 
  if (elecCutHadrEnergy( elec  )!=1){
    return 0; 
  } 
  
  // --- Remove all track conditions on the electron 
  //check for tracks (either bpc or cjc (if in accaptance)): 
  //if(elecCutTrack(elec )!=1) { 
  //  return 0; 
  //}
  
  return 1; 
}


/****************************************************************
 * goodElec_noSpatialCuts 
 *
 * 1. cut on hadronic energy fraction
 * 2. track poiting to cluster 
 */
int elecCut::goodElec_noSpatialCuts( H1PartEm const *elec) { 
 
  //elec whether h1oo finder found it: 
  if(!elec->IsScatElec() ) { 
    return 0; 
  } 

  //make cut on hadronic energy part 
  if (elecCutHadrEnergy( elec  )!=1){
    return 0; 
  } 

  // --- Remove all track conditions on the electron 
  //check for tracks (either bpc or cjc (if in accaptance)): 
  //if(elecCutTrack(elec )!=1) { 
  //  return 0; 
  //}
  return 1; 
}
/****************************************************************
 * goodElec_noSpatialCuts 
 */
int elecCut::goodElec_noSpatialCuts( H1PartEm const *elec, Int_t RunNumber_) { 
  //check if cuts are still ok: 
  initElecCut(RunNumber_); 
  return goodElec_noSpatialCuts(elec); 
}


 /****************************************************************
 * goodElec_NoTrackCond -- all electron cuts exect looking 
 * for a track
 */
int elecCut::goodElec_TrackCond( H1PartEm const *elec) { 
 
  //elec whether h1oo finder found it: 
  if(!elec->IsScatElec() ) { 
    return 0; 
  } 

  //make spatial cuts: 
  if( elecCutSpatial( elec->GetXClus(),  elec->GetYClus() )!=1){ 
    return 0;
  } 
  
  //make cut on hadronic energy part 
  if (elecCutHadrEnergy( elec  )!=1){
    return 0; 
  } 
  
  //check for tracks (either bpc or cjc (if in accaptance)): 
  if(elecCutTrack(elec )!=1) { 
    return 0; 
  }
  
  return 1; 
}

/****************************************************************
 * goodElec_TrackCond -- all electron cuts exect looking 
 * for a track with previous check if cuts are still ok for the
 * given runnumber. 
 */
int elecCut::goodElec_TrackCond( H1PartEm const *elec, Int_t RunNumber_){ 
  initElecCut(RunNumber_); 
  return goodElec_TrackCond(elec); 
} 


/****************************************************************
 * name: elecCutSpatial 
 * parameters:  x and y Position of the Cluster in SpaCal of the 
 *              scattered Electron 
 */
int elecCut::elecCutSpatial(Float_t x, Float_t y ){ 
  
  //Dead Cells: 
  for(int i=0; i<nDeadCells;i++) { 
    if ( x>= posDeadCellsX[i]-DeadCellCut &&   x <= posDeadCellsX[i]+DeadCellCut  &&
	 y>= posDeadCellsY[i]-DeadCellCut &&   y <= posDeadCellsY[i]+DeadCellCut)
      return 0;
  } 

  //Dead Trigger Cells: 
  for(int i=0; i<nDeadTrigCells;i++) { 
    if ( x>= posDeadTrigCellsX[i]-DeadTriggCellCut &&   x <= posDeadTrigCellsX[i]+DeadTriggCellCut  && 
	 y>= posDeadTrigCellsY[i]-DeadTriggCellCut &&   y <= posDeadTrigCellsY[i]+DeadTriggCellCut  )
      return 0;
  } 
  
  //Radius Cut (circle is NOT located in the middle): 
  //Float_t rad =TMath::Sqrt( (x+2.025)*(x+2.025) + y*y ); 
  //if(rad<=RadCut) {
  //  return 0; 
  // }

  //new implementation: 
  Float_t rad; 
  for(int j=0;j<nCircularCuts;j++) { 
    rad =TMath::Sqrt( (x-circCutsXzero[j])*(x-circCutsXzero[j]) + (y-circCutsYzero[j])*(y-circCutsYzero[j]) );
    if(rad<=circCutsRadius[j]) {
      return 0; 
    }
  }

  return 1; 
} 

/****************************************************************
 * elecCutHadrEnergy -- 
 *
 */
int elecCut::elecCutHadrEnergy( H1PartEm const *elec ) { 
  H1FloatPtr ElecE("ElecE");              // energy of scattered electron (acc. to e-finder)
  Float_t EnHadSpac; 

  //the two fraction: 
  Float_t EnFracEmCalo; 
  Float_t EnFracHadCalo; 

  EnFracEmCalo = elec->GetEaem() ;  //energy fraction in electromagnetic part of calo 
  EnHadSpac = elec->GetEnHadSpac();  
  EnFracHadCalo = EnHadSpac / *ElecE;  

  //cut on either fraction (check wheter both are consistent!) 
  if(EnFracHadCalo > MAX_ENERGY_FRAC_HAD) return 0; 

  //energy weighted cluster radius (already made by h1oo finder, just to be able to strengenth cut)
  if(elec->GetEcra() > MAX_ECRA) return 0; 


  return 1; 
} 

/****************************************************************
 * elecCutTrack -- idea is look for any track or hit in the bpc. 
 *                 or in the CJC, depending which acceptance fits
 *                 better. If cluster is neither acceptance range
 *                 accept the event. 
 */
int elecCut::elecCutTrack( H1PartEm const *elec) { 
  Float_t distanceToBPC; 

  //look for any corresponding track: 
  
  //check if in CJC Acceptance: 
  if(*ElecTheta<2.7925) { 
    if(elec->GetTrType()>0) { 
      if(dbg>0) printf("elecCutTrack: Elektron Candidate has corresponding track of type: %d\n", elec->GetTrType()); 
      return 1;   
    }
    else return 0; 
  }

  //check wether we are in the BPC acceptance:
  Double_t x = elec->GetXClus(); 
  Double_t y = elec->GetYClus(); 
  Double_t xshift = 2.5; 
  Double_t r = (TMath::Sqrt( (x+xshift)*(x+xshift) + y*y)); 
  if( (TMath::Abs(y) > 5.5) && r > 19) { 
  
    //in case no track look for bpc hits: 
    distanceToBPC =  elec->GetBcRsp();  // get distance to BPC extrapolation
    if(distanceToBPC < MAX_BPC_DIST) { 
      
      if(dbg>0) printf("elecCutTrack: Elektron has bpc-track at distance %f\n", distanceToBPC); 
      return 1; 
    }
    return 0; 
  }


  if(dbg>0) printf("elecCutTrack: Cluster neither in BPC nor in CJC acceptance.\n");
  return 1; // in neither acceptance -> cannot do something so simply hope that it is no photon
}

/****************************************************************
 * elecCutTrack -- calls above function but checks wether the cuts 
 * are still ok for this event. 
 */
int elecCut::elecCutTrack( H1PartEm const *elec , Int_t RunNumber_) { 
  initElecCut(RunNumber_); 
  return elecCutTrack(elec); 
}


/****************************************************************
 * isInDeadCellArray -- check whether the given position in  
 * the cell which is at the dead cell array at the given position. 
 * Returns 1 if this is the case.  
 */
int elecCut::isInDeadCellArray(H1PartEm const *elec, Int_t i ) { 

    //check range: 
    if(i>=nDeadCells || i<0) { 
	printf("elecCut::isInDeadCellArray: Error request position in dead-cell array is out of bounds\n"); 
	return 0;
    } 
    Float_t x=elec->GetXClus();
    Float_t y=elec->GetYClus(); 
    if(              //cell dimension like for dead trigger cells
	x>= posDeadCellsX[i]- DeadTriggCellCut&&   x <= posDeadCellsX[i]+ DeadTriggCellCut && 
	y>= posDeadCellsY[i]- DeadTriggCellCut&&   y <= posDeadCellsY[i]+ DeadTriggCellCut)
      return 1;
    return 0; 
} 
 

/**********************************************************************
 * DrawCuts
 */ 
void elecCut::DrawCuts( Int_t RunNumber_) { 
  //----- 
  // Draw the cuts into the current pad.
  // Blue lines are dead Trigger Cell Cuts and 
  // red cut are dead cell cut. 
  // -----

  initElecCut(RunNumber_); 
  
  TBox *CellArr[nDeadCells]; 
  TBox *TrigCellArr[nDeadTrigCells]; 
 
 
  //Dead Cells: 
  for(int i=0; i<nDeadCells;i++) { 
    if(dbg) printf("Drawing %d-deadCell(Nmb %d) at: %f %f\n",i,deadCells[i], posDeadCellsX[i], posDeadCellsY[i] ); 

   CellArr[i] = new TBox(posDeadCellsX[i]-DeadCellCut, posDeadCellsY[i]-DeadCellCut,  posDeadCellsX[i]+DeadCellCut,  posDeadCellsY[i]+DeadCellCut); 
    CellArr[i]->SetFillStyle(0); 
    CellArr[i]->SetLineColor(kRed); 
    CellArr[i]->SetLineWidth(2); 
    CellArr[i]->Draw("same"); 
  } 
  
  if(dbg) printf("\n"); 
  //Dead Trigger Cells: 
  for(int i=0; i<nDeadTrigCells;i++) { 
    if(dbg) printf("Drawing %d-deadTriggCell (Nmb %d) at: %f %f\n",i,deadTrigCells[i], posDeadTrigCellsX[i], posDeadTrigCellsY[i] ); 
    TrigCellArr[i] =  new TBox(posDeadTrigCellsX[i]-DeadTriggCellCut, posDeadTrigCellsY[i]-DeadTriggCellCut, posDeadTrigCellsX[i]+DeadTriggCellCut, posDeadTrigCellsY[i]+DeadTriggCellCut); 
      TrigCellArr[i]->SetFillStyle(0); 
    TrigCellArr[i]->SetLineColor(kBlue); 
    TrigCellArr[i]->SetLineWidth(2); 
    TrigCellArr[i]->Draw("same"); 
    
  } 
  

  //Radius Cut (circle is NOT located in the middle): 
  //TEllipse *ell = new TEllipse(-2.025,0, RadCut,RadCut);
  //ell->SetLineColor(kMagenta); ell->SetLineWidth(2);
  //ell->Draw("same");

  //new implementation: 
  if(dbg) printf("Drawing %d circular cuts\n", nCircularCuts); 
  TEllipse *ellArr[nCircularCuts];
  for(int j=0;j<nCircularCuts;j++) { 
    if(dbg) printf("Drawing ciruclar Cut[%d]: at x=%f y=%f with radius=%f\n", j, circCutsXzero[j], circCutsYzero[j], circCutsRadius[j]); 
    //Radius Cut (circle is NOT located in the middle): 
    ellArr[j] = new TEllipse(circCutsXzero[j], circCutsYzero[j], circCutsRadius[j], circCutsRadius[j]);
    ellArr[j]->SetLineColor(kMagenta); ellArr[j]->SetLineWidth(2);
    ellArr[j]->Draw("same");
  }

} 

