#ifndef _ELECCUT_H_
#define _ELECCUT_H_
////////////////////////////////////////////
// Author      : Andreas Jung
// Date        : Aug. 2010
// Description : Implementation of special cuts for scattered electron
//               in SpaCal
////////////////////////////////////////////
#include "H1Mods/H1PartEm.h" 
#include "H1Geom/H1DBManager.h"
#include "H1Pointers/H1FloatPtr.h"

//define of cut thresholds: 
#define MAX_ENERGY_FRAC_HAD    0.15 
#define MAX_ECRA               4     //energy weihted cluster radius 
#define MAX_BPC_DIST           6   // distance to BPC-Hit in cm, to be checked!!!  

//inner radius cut: 
#define RadCut 12  

#define DeadCellCut        3.025  
#define DeadTriggCellCut   2.035
#define MAX_DEAD_CELLS 50 
#define MAX_DEAD_TRIG_CELLS 50
#define MAX_CIRC_CUTS 10 

class elecCut : public TObject { 
private:
  int dbg; 


  Int_t deadCells[MAX_DEAD_CELLS];
  Int_t deadTrigCells[MAX_DEAD_TRIG_CELLS];
  Double_t posDeadCellsX[MAX_DEAD_CELLS]; 
  Double_t posDeadCellsY[MAX_DEAD_CELLS]; 
  Double_t posDeadTrigCellsX[MAX_DEAD_TRIG_CELLS]; 
  Double_t posDeadTrigCellsY[MAX_DEAD_TRIG_CELLS];

  Int_t nCircularCuts; 
  Double_t circCutsRadius[MAX_CIRC_CUTS];
  Double_t circCutsXzero[MAX_CIRC_CUTS]; 
  Double_t circCutsYzero[MAX_CIRC_CUTS]; 

  //initialise the cuts needed for the different run periods: 
  int initElecCut(int RunNumber_) ;
  void clearDeadCellPositions(); 
  void fillDeadCellPositions(); 

  //make the actual cuts: 
  int goodElec_TrackCond( H1PartEm const *elec);
  int elecCutSpatial(Float_t x, Float_t y ); 
  int elecCutHadrEnergy( H1PartEm const *elec ); 
  int elecCutTrack( H1PartEm const *elec); 
  int goodElec_noSpatialCuts( H1PartEm const *elec); 

  Int_t RunNumber; 
  Int_t RunRange; 
  Int_t nDeadCells; 
  Int_t nDeadTrigCells; 
 
public: 
  elecCut(int RunNumber_); 
  elecCut(){}; 
  ~elecCut(){}; 

  void initCuts2004();
  void initCuts2005(); 
  void initCuts2006();
  void initCuts2007();

  //  int initElecCut(int RunNumber); 
  int goodElec( H1PartEm const *elec , Int_t RunNumber_);
  int goodElec_TrackCond( H1PartEm const *elec, Int_t RunNumber_); 
  int goodElec_noSpatialCuts( H1PartEm const *elec, Int_t RunNumber_);
  int elecCutTrack( H1PartEm const *elec , Int_t RunNumber_); 
  void DrawCuts( Int_t RunNumber_);
  int isInDeadCellArray(H1PartEm const *elec, Int_t i ); 

  //ClassDef(elecCut,1) 

}; 
#endif
