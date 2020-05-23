// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch, based on main01.cc.
// This program illustrates how HepMC can be interfaced to Pythia8.
// It studies the charged multiplicity distribution at the LHC.
// HepMC events are output to the hepmcout41.dat file.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include <math.h>
#include <algorithm>
#include <list>
#include "fastjet/PseudoJet.hh"
using namespace HepMC;


int Mother(HepMC::GenParticle* p, int n){
     if ( (p)->production_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->production_vertex()->particles_begin(HepMC::parents);
    HepMC::GenVertex::particle_iterator up = (p)->production_vertex()->particles_end(HepMC::parents);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[10]={0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
         x[i]= (*mother)->barcode();
        i++;
        }
         return x[n];
         delete x;
    }
//     else{return -1;}   
}


int Mother1(HepMC::GenParticle* p){
     if ( (p)->production_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->production_vertex()->particles_begin(HepMC::parents);
    HepMC::GenVertex::particle_iterator up = (p)->production_vertex()->particles_end(HepMC::parents);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[2]={0,0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
         x[i]= (*mother)->barcode();
        i++;
        }
         return x[0];
         delete x;
    }
//     else{return -1;}   
}

int Mother2(HepMC::GenParticle* p){
     if ( (p)->production_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->production_vertex()->particles_begin(HepMC::parents);
    HepMC::GenVertex::particle_iterator up = (p)->production_vertex()->particles_end(HepMC::parents);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[2]={0,0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
         x[i]= (*mother)->barcode();
        i++;
        }
         return x[1];
         delete x;
    }
//     else{return -1;}   
}

int Daughter(HepMC::GenParticle* p, int n){
     if ( (p)->end_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
    HepMC::GenVertex::particle_iterator up = (p)->end_vertex()->particles_end(HepMC::children);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[10]={0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
         x[i]= (*mother)->barcode();
        i++;
        }
         return x[n];
         delete x;
    }
//     else{return -1;}   
}

//if ( end_vertex() && end_vertex()->barcode()!=0 )
int Daughter1(HepMC::GenParticle* p){
     if ( (p)->end_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
    HepMC::GenVertex::particle_iterator up = (p)->end_vertex()->particles_end(HepMC::children);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[2]={0,0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator daughter=lp; daughter!=up;++daughter) { // particle loop for each event 
         x[i]= (*daughter)->barcode();
        i++;
        }
         return x[0];
         delete x;
    }
//     else{return -1;}   
}

int Daughter2(HepMC::GenParticle* p){
     if ( (p)->end_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
    HepMC::GenVertex::particle_iterator up = (p)->end_vertex()->particles_end(HepMC::children);
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();
    int x[2]={0,0};
    int i = 0;
    for ( HepMC::GenVertex::particle_iterator daughter=lp; daughter!=up;++daughter) { // particle loop for each event 
         x[i]= (*daughter)->barcode();
        i++;
        }
         return x[1];
         delete x;
    }
//     else{return -1;}   
}


inline bool IsParton(HepMC::GenParticle* p){
    return (abs(p->pdg_id())<9||p->pdg_id()==21 );
}

inline bool IsFinal(HepMC::GenParticle* p){
    return (!p->end_vertex()&&p->status()==1 );
}

inline bool IsNeutrino(HepMC::GenParticle* p){
    return ((abs(p->pdg_id())- 12)*(abs(p->pdg_id())- 14)*(abs(p->pdg_id())- 16) ==0 );
}

inline bool IsHadron(HepMC::GenParticle* p){
    int idSave = abs(p->pdg_id());
    if (idSave <= 100 || (idSave >= 1000000 && idSave <= 9000000)|| idSave >= 9900000){ return false;}
    if (idSave == 130 || idSave == 310){ return true;}
    if (idSave%10 == 0 || (idSave/10)%10 == 0 || (idSave/100)%10 == 0){return false;}
  return true;
}

inline bool IsPhoton(HepMC::GenParticle* p){
    return (p->pdg_id()== 22 );
}



int main() {

int nEvent=2000;    

  // specify an input file
// HepMC::IO_GenEvent ascii_in("/home/ja2006203966/event/HerwigLHC.hepmc",std::ios::in);
// HepMC::IO_GenEvent ascii_in("/home/ja2006203966/event/Sherpa/sherpahep.hepmc2g",std::ios::in);
HepMC::IO_GenEvent ascii_in("/home/ja2006203966/event/test.hepmc",std::ios::in);

// get the first event
HepMC::GenEvent* evt = ascii_in.read_next_event();
HepMC::GenEvent::particle_iterator pit;
	int icount=0;
	int num_good_events=0;
// loop until we run out of events
// for ( HepMC::GenEvent::particle_const_iterator p = HepMC::particles_begin();
// 	           p != HepMC::particles_end(); ++p ) {
// 	        (*p)->print();
// 	     }

// just read parents's id
int iEvent=0;
int i=0;
// std::cout<<evt;
while(evt){
     if(i==10){break;}
     std::cout<<i<<"-th event\n";
    
 for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();p != evt->particles_end(); ++p ) {
// 	    (*p)->print();
//     std::cout<<"no.\t"<<(*p)->barcode() <<"\tp4\t"<<(*p)->momentum().perp()<<"\n";
    
// HepMC::GenEvent mother = (*p)->parent_event() ;
    
//      std::cout<<"particle number\t"<<(*p)->barcode()<<"\tparticle id\t"<<(*p)->pdg_id()<<"\n";
//     std::cout<<"particle number\t"<<(*p)->barcode()<<"\t first particle id\t"<<evt->barcode_to_particle(5)->pdg_id()<<"\n";
  if (IsNeutrino(*p)){
    std::cout<<"particle number\t"<<(*p)->barcode()<<" particle id\t"<<(*p)->pdg_id()<<"\t mother1 num\t"<<Mother1(*p)<<"\t mother2 num\t"<<Mother2(*p)<<
        "\t daughter1 num\t"<<Daughter1(*p)<<"\t daughterr2 num\t"<<Daughter2(*p)<<"\n";
      
    if (IsParton(*p)){std::cout<<"IsParton\n";}
    
    if (IsFinal(*p)){std::cout<<"IsFinal\n";}
    
    if (IsHadron(*p)){std::cout<<"IsHadron\n";}
      std::cout<<"mother1 id "<<evt->barcode_to_particle(Mother1(*p))->pdg_id()<<"\n";
      
    
//     int test[10] = {0};
//     for (int i: test){
//         std::cout<<"test"<<i<<"\n";
//     }
//     for (int i; i<10;i++){
//         std::cout<<Mother(p,i)<<"\n";
//     }
      
//     std::cout << "\t";
//     HepMC::GenEvent::particle_iterator p = evt->particles_begin();
    
    if ( (*p)->production_vertex() ) {
    HepMC::GenVertex::particle_iterator lp = (*p)->production_vertex()->particles_begin(HepMC::parents);
    HepMC::GenVertex::particle_iterator up = (*p)->production_vertex()->particles_end(HepMC::parents);
    HepMC::GenVertex::particle_iterator mother=lp;
//     std::cout << "\t";
//     std::cout<<(*mother)->pdg_id();    
    for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) {
        std::cout << "\t";
        std::cout<<"mother num \t"<<(*mother)->barcode()<<"\t";
        std::cout<<"mother id "<<(*mother)->pdg_id()<<"\n";
                    }
        
	}
  }
    if (iEvent==nEvent){break;}
 iEvent++;   
 }
i++;
delete evt;
evt = ascii_in.read_next_event();
}


  return 0;
}
