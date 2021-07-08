// (c) MPI 2020
#include "H2020HistManager.h"
 
// C++ includes
#include <fstream>
#include <set>

// Root includes
 
// H1 includes

using namespace std;


// _______________________________________________________ //
//! Constructor
H2020HistManager::H2020HistManager(const string& HMname, const string& dirname) : fHMname(HMname) {
   TH1::AddDirectory(false);
   fDirname = dirname;
   if ( fDirname == "" ) fDirname = fHMname;
}

// _______________________________________________________ //E
H2020HistManager::~H2020HistManager() {

}

// _______________________________________________________ //
//! Write histograms to gDirectory
void H2020HistManager::Write() {
   cout<<"H2020HistManager::Write. Writing histograms into "<<gDirectory->GetPath()<<endl;

   // sort by name first
   map<string,TH1*> hmap;
   for ( auto [k1, m1] : fHistmap ) {
      for ( auto [k2,th] : m1 ) {
         hmap[th->GetName()] = th;
      }
   }
   for ( auto [k1, th] : hmap )  th->Write();
}