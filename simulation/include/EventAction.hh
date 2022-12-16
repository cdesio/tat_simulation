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
#include <vector>
#include "G4KDTree.hh"
#include "DetectorConstruction.hh"

class DetectorConstruction; 

class EventAction : public G4UserEventAction
{
public:
  EventAction(DetectorConstruction* pDetector);
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event *);
  virtual void EndOfEventAction(const G4Event *);

  G4double GetEdep() { return fEdep; };
  void SetEdep(G4double pVal) { fEdep = pVal; };
  void AddEdep(G4double pVal) { fEdep += pVal; };

  // Record the track length though the chomatin fibre segment
  void SetStartTrackKE(G4double pVal) { fTrackStartKE = pVal; };
  void SetStartTrackFound() { fTrackStartFound = true; };
  G4bool GetStartTrackFound() { return fTrackStartFound; };

  void SetEndTrackKE(G4double pVal) { fTrackEndKE = pVal; };

  G4ThreeVector GetStartTrackPos() { return fTrackStartPos; };
  void SetStartTrackPos(G4ThreeVector pVal) { fTrackStartPos = pVal; };

  G4ThreeVector GetEndTrackPos() { return fTrackEndPos; };
  void SetEndTrackPos(G4ThreeVector pVal) { fTrackEndPos = pVal; };

  G4bool GetStoppedInBox() { return fTrackStoppedBox; };
  void SetStoppedInBox() { fTrackStoppedBox = true; };

  void AddPathLength(G4double val) { fpathLengthTotal += val; }
  G4double Getdkill(){return dkill;}

  std::vector<G4double> fTrackMeanKE;
  G4KDTree *fPositions0Event;
  G4KDTree *fPositions1Event;
  G4KDTree *fPositionsBase0Event;
  G4KDTree *fPositionsBase1Event;

private:
  G4double fEdep;
  G4double fTrackStartKE;
  G4double fTrackEndKE;
  G4bool fTrackStartFound;
  G4ThreeVector fTrackStartPos;
  G4ThreeVector fTrackEndPos;
  G4bool fTrackStoppedBox;
  G4double fpathLengthTotal{0};
  G4double dkill;

};

#endif