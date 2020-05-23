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
#include <iostream>
#include <cmath>
#include <fstream>

// head file for jet clustering
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/DistanceMeasure.hh"
#include "fastjet/contrib/QCDAwarePlugin.hh"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include <math.h>
#include <algorithm>
#include <list>

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



double DeltaPhi(double phi1, double phi2){
  if (abs(phi1 - phi2) > M_PI) return 2 * M_PI - abs(phi1 - phi2);
  else return abs(phi1 - phi2);
}

double DeltaR(double phi1, double phi2, double eta1, double eta2){
  return sqrt(pow(DeltaPhi(phi1, phi2), 2) + pow(eta1 - eta2, 2));
}


int main(int argc, char *argv[]) {

  // specify an input file
HepMC::IO_GenEvent ascii_in("/home/ja2006203966/event/test.hepmc",std::ios::in);
// get the first event
HepMC::GenEvent* evt = ascii_in.read_next_event();
HepMC::GenEvent::particle_iterator pit;
using namespace std;

// loop until we run out of events
// for ( HepMC::GenEvent::particle_const_iterator p = HepMC::particles_begin();
// 	           p != HepMC::particles_end(); ++p ) {
// 	        (*p)->print();
// 	     }

// just read parents's id
 ofstream myfile;
   myfile.open ("/home/ja2006203966/event/HepmctoML.txt"); 
//============================================Fastjet setting================================

 // Fastjet analysis - select algorithm and parameters.
  double Rparam = 0.4;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition *jetDef = NULL;
  jetDef = new fastjet::JetDefinition( fastjet::antikt_algorithm, Rparam,recombScheme, strategy);
  int Maxpt = 550; // maxima jet pt
  int Minpt = 500; // minima jet pt
// QCD - aware Fastjet analysis
  double Ghostparam = 1e-20;
  fastjet::contrib::QCDAwarePlugin::AntiKtMeasure *akt = new fastjet::contrib::QCDAwarePlugin::AntiKtMeasure(Rparam);
  fastjet::contrib::QCDAwarePlugin::QCDAwarePlugin *qcdawareakt = new fastjet::contrib::QCDAwarePlugin::QCDAwarePlugin(akt);

  // Fastjet input.
  std::vector <fastjet::PseudoJet> fjInputs, QCDfjInputs;

  // set parameter for jet splitting
  double dR = 0.2;
  int nEvent = 100;

  double MassMax = 10000.;
  double pTMax = 5000.;
  int iEvent = 0;

 // hard coding to find 2 partons before showering
 
    fjInputs.clear();

//====================================HepMC2 loop all event====================================================
while(evt){ // loop all event
    
    std::cout<<iEvent<<"-th event\n";
//======================================QCD-aware=======================================

  // Fastjet analysis - select algorithm and parameters.
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition( fastjet::antikt_algorithm, Rparam,
           recombScheme, strategy);

  // QCD - aware Fastjet analysis
  double Ghostparam = 1e-20;
  fastjet::contrib::QCDAwarePlugin::AntiKtMeasure *akt = new fastjet::contrib::QCDAwarePlugin::AntiKtMeasure(Rparam);
  fastjet::contrib::QCDAwarePlugin::QCDAwarePlugin *qcdawareakt = new fastjet::contrib::QCDAwarePlugin::QCDAwarePlugin(akt);

  // Fastjet input.
  std::vector <fastjet::PseudoJet> fjInputs, QCDfjInputs;

  // set parameter for jet splitting
  double dR = 0.2;

  int nEvent = 100;
  double MassMax = 10000.;
  double pTMax = 5000.;

  // output file with tried event number & cross section
 
    ofstream information;
  if (argc > 3) information.open(argv[3], ios::out);
  else information.open("Info.txt", ios::out);

  // settings for using ROOT histograms
 

 
    // hard coding to find 2 partons before showering
    
    fjInputs.clear();
//=========================Particle(expect neutrino)  loop  for particle jet===========================
    //HepMC2 version
    for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();p != evt->particles_end(); ++p ) {
        // if you want Final particle : 
        
//==================================================================
      // No neutrinos or DM.
      // 12 for electron neutrino, 14 for muon neutrino, 16 for tauon neutrino, 52 for dark matter with spin 1 / 2
      // Pdgid can be accessed in https://twiki.cern.ch/twiki/bin/view/Main/PdgId
      if (!IsFinal(*p)){continue;}
      if (IsNeutrino(*p)){continue;}

      // Only |eta| < 3.6.
      //      if (abs(pythia.event[i].eta()) > 3.6) continue;
      // still don't know how to cut.
      // put every particle into list?

      // Store as input to Fastjet.
      fjInputs.push_back(fastjet::PseudoJet( (*p)->momentum().px(),(*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e() ) );
    }
//========================================end of Particle(expect neutrino) loop===================================================
    

    // Check that event contains analyzable particles.
    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm.
    std::vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    // need at least 2 jets to finish leading jet and sub-leading jet analysis
    if (sortedJets.size() < 2) {
      cout << "No enough jets found in event " << iEvent << endl;
      continue;
    }
/* ----------------pythia version -> hepmc2  (Vec4::Vec4(double x = 0., double y = 0., double z = 0., double t = 0.) to  	FourVector (double xin, double yin, double zin, double tin=0) )-------------*/
   /* Vec4 pJ1(sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].e());
    Vec4 pJ2(sortedJets[1].px(), sortedJets[1].py(), sortedJets[1].pz(), sortedJets[1].e());
    Vec4 pC = pJ1 + pJ2; */
    HepMC::FourVector pJ1(sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].e());
    HepMC::FourVector pJ2(sortedJets[1].px(), sortedJets[1].py(), sortedJets[1].pz(), sortedJets[1].e());
    HepMC::FourVector pC(pJ1.px()+pJ2.px(), pJ1.py()+pJ2.py(),pJ1.pz()+pJ2.pz(), pJ1.pz()+pJ2.pz()) ;
//----------------------------------------------------
    double Ystar = (sortedJets[0].rap() - sortedJets[1].rap()) / 2;

    // cut events with too soft leading jets and sub leading jets
    //------------ .pT() -> .per()
    if ((pJ1.perp() < Maxpt) || (pJ2.perp() < Minpt)){
      cout << "Jets too soft in event " << iEvent << endl; 
      continue;
    }

    // cut the events where jet eta is too large
    if ((abs(pJ1.eta()) >= 2) || (abs(pJ2.eta()) >= 2)){
      cout << "Jets with too large eta in event " << iEvent << endl;
      continue;
    }

    // Qcdaware jet clustering

    QCDfjInputs.clear();
//=================================== particle (pick final parton) loop particle jet====================================================
    // Particle loop to pick final partons
    for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();p != evt->particles_end(); ++p ) {
      // partons only
      if (!IsParton(*p)){ continue;}
      // no parton child
     // evt->barcode_to_particle(Daughter1(p))
      if (IsParton(evt->barcode_to_particle(Daughter1(*p))) || IsParton(evt->barcode_to_particle(Daughter2(*p)))){continue;}
      // reject if the parton is from hadron or tau decay
      bool fromHorT = 0;
      int mother ;
        
      for (int test; test<10; test++){//check all mother's status =2 and if have least one is hadron
        mother = Mother(*p,test);
        if (mother=0){continue;}
          fromHorT = fromHorT || (IsHadron(evt->barcode_to_particle(mother)) && evt->barcode_to_particle(mother)->status() == 2);
          fromHorT = fromHorT || ((evt->barcode_to_particle(mother)->pdg_id() == 15) && ( evt->barcode_to_particle(mother)->status() == 2));
        }
        
      if (fromHorT){continue;} // if above(at least on mother is hadron and status =2) continue

      // put the selected parton in fastjet input vector
      fastjet::PseudoJet FinalParton((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
      FinalParton.set_user_index((*p)->pdg_id());

      QCDfjInputs.push_back(FinalParton);
    }
//====================================end of particle (pick final parton) loop=================================================================

    // if no final partons
    if (fjInputs.size() == 0){
      cout << "Error: no final partons in event " << iEvent << endl;
      continue;
    }

    // QCD aware jet clustering
    std::vector <fastjet::PseudoJet> QCDSortedJets;

    fastjet::ClusterSequence QCDclustSeq(QCDfjInputs, qcdawareakt);
    QCDSortedJets = sorted_by_pt(QCDclustSeq.inclusive_jets(20.0));

    // no enough parton jets
    if (QCDSortedJets.size() < 2){
      cout << "No enough parton jets found in event" << iEvent << endl;
      continue;
    }

    // add "ghost" parton jet to the pseudojet list
    // ghostify the parton jet
    std::vector <fastjet::PseudoJet> refjInput;
    refjInput.clear();
    // default user index of pseudo jet is -1, the same as anti-down quark
    // so set the user index of particle jet as 0
    for (fastjet::PseudoJet particlejet: sortedJets){
      particlejet.set_user_index(0);
      refjInput.push_back(particlejet);
    }
    for (fastjet::PseudoJet partonjet: QCDSortedJets){
      // ghostify
      //      cout << "Partonjet:"
      //	   << partonjet.eta() << " "
      //	   << partonjet.phi() << " ";
      partonjet.reset_momentum(partonjet.px() * Ghostparam, partonjet.py() * Ghostparam, partonjet.pz() * Ghostparam, partonjet.e() * Ghostparam);
      //      cout << partonjet.pt() / Ghostparam << " "
      //	   << partonjet.e() / Ghostparam << " "
      //	   << partonjet.user_index() << " "
      //	   << endl;
      refjInput.push_back(partonjet);
    }

    // reclustering
    fastjet::ClusterSequence reClustSeq(refjInput, *jetDef);
    std::vector <fastjet::PseudoJet> reSortedJets;
    reSortedJets = sorted_by_pt(reClustSeq.inclusive_jets(20.0));

    // labelling
    // pick the closest parton jets to label
    // if the distance of all the parton constituents are larger than dR, the particle jet remain unlabelled
    int EventType = -1;
    double Rmin0 = dR, Rmin1 = dR;
    int label0 = 0, label1 = 0;
    for (fastjet::PseudoJet con: reSortedJets[0].constituents())
      if (con.user_index() != 0)
	if (DeltaR(reSortedJets[0].phi(), con.phi(), reSortedJets[0].eta(), con.eta()) < Rmin0){
	  Rmin0 = DeltaR(reSortedJets[0].phi(), con.phi(), reSortedJets[0].eta(), con.eta());
	  label0 = con.user_index();
	}
    for (fastjet::PseudoJet con: reSortedJets[1].constituents())
      if (con.user_index() != 0)
	if (DeltaR(reSortedJets[1].phi(), con.phi(), reSortedJets[1].eta(), con.eta()) < Rmin1){
	  Rmin1 = DeltaR(reSortedJets[1].phi(), con.phi(), reSortedJets[1].eta(), con.eta());
	  label1 = con.user_index();
	}

    // Delta R are all larger than dR
    if ((label0 == 0) || (label1 == 0)){
      
      continue;
    }
    label0 = abs(label0);
    label1 = abs(label1);
//==================================================record data===============================================
//      if ((label0 <= 8) && (label1 <= 8)){
//       myfile<<"\n"<<sortedJets[i].e()<<"\t"<<sortedJets[i].pt()<<"\t"<<sortedJets[i].eta()<<"\t"<<sortedJets[i].phi()<<"\t"<<constituents.size()<<"\n"<<endl; 
// 			for (int j=0;j<int(constituents.size()); ++j){  
// myfile<<constituents[j].e()<<"\t"<<constituents[j].pt()<<"\t"<<constituents[j].eta()<<"\t"<<constituents[j].phi()<<endl;
//     }
//      }
//     else if ((label0 == 21) && (label1 == 21)){
//       HMassCgg->Fill(pC.mCalc());
//       HDijetType->Fill(1);
//     } 
//    else if (((label0 == 21) && (label1 <= 8)) ||((label0 <= 8) && (label1 == 21))){
//       HMassCqg->Fill(pC.mCalc());
//       HDijetType->Fill(2);
//     } 
//     else {
//       HMassCun->Fill(pC.mCalc());
//       HDijetType->Fill(-1);
//     }
  
 
  delete evt;
evt = ascii_in.read_next_event();
iEvent++;
}

myfile.close();
  return 0;
}
