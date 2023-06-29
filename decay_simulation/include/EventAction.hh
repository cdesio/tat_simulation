//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"


class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  G4double GetEdep() {return fEdep;};
  void SetEdep(G4double pVal) {fEdep=pVal;};
  void AddEdep(G4double pVal) {fEdep+=pVal;};


  void AddSecondary(){numSecondary+=1;}
  G4int GetNumSecondaries(){return numSecondary;}

  // Record the path length though the chomatin fibre segment
  void AddPathLength(G4double val){fpathLengthTotal += val;}
  std::map<G4int, G4int> parentParticle;

private:
  G4double  fEdep;
  G4int numSecondary;
  G4double fpathLengthTotal{0};
  std::map<G4int, G4String> particleMapRev{
      {1, "alpha"},
      {2, "gamma"},
      {3, "e-"},
      {4, "nu_e"},
      {5, "At211"},
      {6, "Po211"},
      {7, "Bi207"},
      {8, "Pb207"},
      {9, "e+"}};
};

#endif