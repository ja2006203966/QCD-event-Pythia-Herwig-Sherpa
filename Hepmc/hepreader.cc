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
         if(n<=(p->production_vertex()->particles_in_size())){
             HepMC::GenVertex::particle_iterator lp = (p)->production_vertex()->particles_begin(HepMC::parents);
             for(int i=1; i<n; ++i){ lp = ++lp;}
             return (*lp)->barcode();
         }
         if(n>(p->production_vertex()->particles_in_size())){return 0;}
   
//     HepMC::GenVertex::particle_iterator up = (p)->production_vertex()->particles_end(HepMC::parents);
    
//     int x[10]={0};
//     int i = 0;
//     for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
//          x[i]= (*mother)->barcode();
//         i=i+1;
//         }
//          return x[n];
//          delete x;
    }
    else{return 0;}   
}


int Mother1(HepMC::GenParticle* p){
     if ( p->production_vertex() ) {
         HepMC::GenVertex::particle_iterator lp = p->production_vertex()->particles_begin(HepMC::parents);
         return (*lp)->barcode();
    }
    else{return 0;}   
}

int Mother2(HepMC::GenParticle* p){
     if ( p->production_vertex() ) {
         if(p->production_vertex()->particles_in_size()>=2){
             HepMC::GenVertex::particle_iterator lp = (p)->production_vertex()->particles_begin(HepMC::parents);
             lp = ++lp;
         return (*lp)->barcode();
         }
         else{return 0;}
    }
    else{return 0;}   
}

int Daughter(HepMC::GenParticle* p, int n){
     if ( p->end_vertex() ) {
         if(n<=(p->end_vertex()->particles_out_size())){
             HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
             for(int i=1; i<n; ++i){ lp = ++lp;}
             return (*lp)->barcode();
         }
         else{return 0;}

    }
    else{return 0;}   
}

// int Daughter(HepMC::GenParticle* p, int n){
//      if ( (p)->end_vertex() ) {
//     HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
//     HepMC::GenVertex::particle_iterator up = (p)->end_vertex()->particles_end(HepMC::children);
// //     std::cout << "\t";
// //     std::cout<<(*mother)->pdg_id();
//     int x[10]={0};
//     int i = 0;
//     for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) { // particle loop for each event 
//          x[i]= (*mother)->barcode();
//         i++;
//         }
//          return x[n];
// //          delete x;
//     }
// //     else{return -1;}   
// }

//if ( end_vertex() && end_vertex()->barcode()!=0 )
int Daughter1(HepMC::GenParticle* p){
     if ( (p)->end_vertex() ) {
         HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
         return (*lp)->barcode();
    }
    else{return 0;}   
}

int Daughter2(HepMC::GenParticle* p){
     if ( (p)->end_vertex() ) {
         if(p->end_vertex()->particles_out_size()>=2){
             HepMC::GenVertex::particle_iterator lp = (p)->end_vertex()->particles_begin(HepMC::children);
             lp = ++lp;
             return (*lp)->barcode();
         }
         else{return 0;}
       
    }
    else{return 0;}   
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
// HepMC::GenEvent* evt = ascii_in.read_next_event();
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
// while(evt){

for (HepMC::GenEvent* evt = ascii_in.read_next_event(); evt; evt=ascii_in.read_next_event()){
    
     if(i==1){break;}
    i=i+1;
     std::cout<<evt->event_number()<<"-th event\t"<<"event size=\t"<<evt->particles_size()<<"\n";
//     std::cout<<"914 id:\t"<<evt->barcode_to_particle(914)->pdg_id()<<"\n";
    int check = 0;
 
 for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();p != evt->particles_end(); ++p ) {
// 	    (*p)->print();
//     std::cout<<"no.\t"<<(*p)->barcode() <<"\tp4\t"<<(*p)->momentum().perp()<<"\n";
    
// HepMC::GenEvent mother = (*p)->parent_event() ;
    
//      std::cout<<"particle number\t"<<(*p)->barcode()<<"\tparticle id\t"<<(*p)->pdg_id()<<"\n";
//     std::cout<<"particle number\t"<<(*p)->barcode()<<"\t first particle id\t"<<evt->barcode_to_particle(5)->pdg_id()<<"\n";
//   if (true){
     
//      (*p)->print();
//      int m2 = Mother2(*p);
    std::cout<<"particle number\t"<<(*p)->barcode()<<" particle id\t"<<(*p)->pdg_id()<<"\t mother1 num\t"<<Mother1(*p)<<"\t mother2 num\t"<<Mother2(*p)<<"\t daughter1 num\t"<<Daughter1(*p)<<"\t daughterr2 num\t"<<Daughter2(*p)<<"\n";
     
     for(int mother=1;mother!=5;mother++){
         std::cout<<mother<<"-th mother =\t"<<Mother(*p,mother)<<"\t";
     }
     std::cout<<"\n";
     
     for(int daughter=1;daughter!=5;daughter++){
         std::cout<<daughter<<"-th daughter =\t"<<Daughter(*p,daughter)<<"\t";
     }
     std::cout<<"\n";
     
    if (IsParton(*p)){std::cout<<"IsParton\n";}
    
    if (IsFinal(*p)){std::cout<<"IsFinal\n";}
    
    if (IsHadron(*p)){std::cout<<"IsHadron\n";} 
     
     
//      if(Mother1(*p)!=0){ std::cout<<"mother1 id "<<evt->barcode_to_particle(Mother1(*p))->pdg_id()<<"\n";}
      
    
//     int test[10] = {0};
//     for (int i: test){
//         std::cout<<"test"<<i<<"\n";
//     }
//     for (int i; i<10;i++){
//         std::cout<<Mother(p,i)<<"\n";
//     }
      
//     std::cout << "\t";
//     HepMC::GenEvent::particle_iterator p = evt->particles_begin();
    
     
     
     
     
     //=================================== particle (pick final parton) loop particle jet====================================================

 
//       // partons only
      
//       if (!IsParton(*p)){check++; continue;}
// //         std::cout<<iEvent<<"-th event\n check core not dump6\n";
//       // no parton child
//      // evt->barcode_to_particle(Daughter1(p))
//         int d1=Daughter1(*p);
//         int d2=Daughter2(*p);
         
//         if(d1!=0){
           
//           if (IsParton(evt->barcode_to_particle(d1)) || IsParton(evt->barcode_to_particle(d2))){check++; continue;}
//          }
         
//       // reject if the parton is from hadron or tau decay
//       bool fromHorT = 0;
//       int mother ;
        
        
//       for (int test=0; test<10; test++){//check all mother's status =2 and if have least one is hadron
//         mother = Mother(*p,test);

//         if (mother==0){continue;}
        
// //           if(IsHadron(evt->barcode_to_particle(mother)) && evt->barcode_to_particle(mother)->status() == 2){fromHorT=1;}
// //           if((evt->barcode_to_particle(mother)->pdg_id() == 15)&& evt->barcode_to_particle(mother)->status() == 2){fromHorT=1;}
          
// //           eq1 = IsHadron(evt->barcode_to_particle(mother));
// //           std::cout<<iEvent<<"-th event\n check core not dump10\n";
// //           eq2 = (evt->barcode_to_particle(mother)->status() == 2);
// //           std::cout<<iEvent<<"-th event\n check core not dump11\n";
// //           eq3 = (evt->barcode_to_particle(mother)->pdg_id() == 15);
// //           std::cout<<iEvent<<"-th event\n check core not dump12\n";       
// //           fromHorT = fromHorT || IsHadron(evt->barcode_to_particle(mother)) && (evt->barcode_to_particle(mother)->status() == 2);
// //           fromHorT = fromHorT || (evt->barcode_to_particle(mother)->pdg_id() == 15) && (evt->barcode_to_particle(mother)->status() == 2);
//         }

  
// //       if (fromHorT){ check++; continue;} // if above(at least on mother is hadron and status =2) continue

//       // put the selected parton in fastjet input vector
//       fastjet::PseudoJet FinalParton((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
//       FinalParton.set_user_index((*p)->pdg_id());
//         std::cout<<iEvent<<"-th event check core not dump2\n";

//====================================end of particle (pick final parton) loop=================================================================
     
     
     
     
     
     
     
     
     
     
     
     
     
     
//     if ( (*p)->production_vertex() ) {
//     HepMC::GenVertex::particle_iterator lp = (*p)->production_vertex()->particles_begin(HepMC::parents);
//     HepMC::GenVertex::particle_iterator up = (*p)->production_vertex()->particles_end(HepMC::parents);
// //     std::cout << "\t";
// //     std::cout<<(*mother)->pdg_id();    
//     for ( HepMC::GenVertex::particle_iterator mother=lp; mother!=up;++mother) {
//         std::cout << "\t";
//         std::cout<<"mother num \t"<<(*mother)->barcode()<<"\t";
//         std::cout<<"mother id "<<(*mother)->pdg_id()<<"\n";
//                     }
        
// 	}
//   }
     
    
//     if (iEvent==nEvent){break;}
 iEvent=iEvent+1;   
 }

// delete evt;
// evt = ascii_in.read_next_event();
}


  return 0;
}
