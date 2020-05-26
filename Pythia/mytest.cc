// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8Plugins/FastJet3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <set>
#include <stdlib.h> 

#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"


using namespace std;
using namespace Pythia8;


void print(std::vector<int> const &input)
{
	for (auto const& i: input) {
		std::cout << i << " ";
	}
}

int main(int argc, char *argv[]) {
  
  int  nEvent = 50000;
  int nListJets = nEvent; //this value should <= nEvent
 
  Pythia pythia;
  
  Event& event = pythia.event;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");//"Random:seed = 0"
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 12 -12");
  pythia.readString("23:onIfMatch = 14 -14");
  pythia.readString("23:onIfMatch = 16 -16");
  pythia.readString("PhaseSpace:pTHatMin = 500.");  //PhaseSpace:pTHatMax   (default = -1.)
  pythia.readString("PhaseSpace:pTHatMax = 550.");
  
  

  pythia.init();
//========================================================================================
  //open file

  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io ("/home/ja2006203966/event/test3.hepmc", std :: ios :: out );
//ofstream myout;
  //myout.open ("myout.txt");//qq for quark ,gg for gluon
//========================================================================================


  // Begin event loop. Generate event. Skip if error.
for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
   
   
    if (!pythia.next()) continue;
     

ToHepMC.fill_next_event( pythia, hepmcevt );
ascii_io << hepmcevt;
delete hepmcevt;



  
} 
 //=================================================================
   //close file
 
 // myout.close();
  return 0;
}
