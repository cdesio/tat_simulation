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
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include <G4Scheduler.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAMolecule.hh"
#include "G4MoleculeTable.hh"
#include "G4OH.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TransportationManager.hh"
#include "CommandLineParser.hh"
// #include "G4LogicalVolumeStore.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4KDTreeResult.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(DetectorConstruction *pDetector)
    : G4UserTimeStepAction(), fpDetector(pDetector)
{

  EventAction *fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
  fPositions0 = fpEventAction->fPositions0Event;
  fPositions1 = fpEventAction->fPositions1Event;
  fPositionsBase0 = fpEventAction->fPositionsBase0Event;
  fPositionsBase1 = fpEventAction->fPositionsBase1Event;
  dkill = ((DetectorConstruction *)fpDetector)->Getdkill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction &other)
    : G4UserTimeStepAction(other)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction &
TimeStepAction::operator=(const TimeStepAction &rhs)
{
  if (this == &rhs)
    return *this;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void TimeStepAction::UserPostTimeStepAction()
{
  G4TrackManyList *allTrackList = G4ITTrackHolder::Instance()->GetMainList();
  G4TrackManyList::iterator it = allTrackList->begin();
  G4TrackManyList::iterator end = allTrackList->end();
  G4Track *trackToKill;

   while (it != end)

  {
    trackToKill = nullptr;

    G4Track *track = *it; // track can be an OH, e_aq, H2, ...
    G4String name = GetMolecule(track)->GetName();

    if ((name != "Deoxyribose^0") && (name != "Cytosine^0") && (name != "Adenine^0") && (name != "Thymine^0") && (name != "Guanine^0"))
    {
      auto localPos = track->GetPosition();

      // Remove all radicals more than dkill from DNA
      G4KDTreeResultHandle result0 = fPositions0->NearestInRange(localPos, dkill);
      G4KDTreeResultHandle result1 = fPositions1->NearestInRange(localPos, dkill);

      if ((result0->GetSize() == 0) && (result1->GetSize() == 0))
      {
        G4ITReactionSet *fReactionSet = G4ITReactionSet::Instance();
        fReactionSet->RemoveReactionSet(track);
        // track->SetTrackStatus(fKillTrackAndSecondaries);
              trackToKill = track;


        // continue;
      }

      G4VTouchable *newTouchable = CreateTouchableAtPoint(localPos);

      G4String volumeName = newTouchable->GetSolid()->GetName();

      if ((volumeName == "histone") || (volumeName == "world") || (volumeName == "waterBox"))

      {
        // Histones act as perfect scavengers for all radiolysis products
        // Do not want to track radiolysis products in the other volumes to save computational time, this leaves volumes only immediately around the bases.
        G4ITReactionSet *fReactionSet = G4ITReactionSet::Instance();
        fReactionSet->RemoveReactionSet(track);
        // track->SetTrackStatus(fKillTrackAndSecondaries);
        // G4ITTrackHolder::Instance()->PushToKill(track);
      trackToKill = track;


      }

      delete newTouchable;

    }
    ++it;
    if (trackToKill != nullptr)
    {
      // need to remove all reactions containing this track for IRT method
      G4ITReactionSet* fReactionSet = G4ITReactionSet::Instance();
      fReactionSet->RemoveReactionSet(trackToKill);
      G4ITTrackHolder::Instance()->PushToKill(trackToKill);
    }
  }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserReactionAction(const G4Track &trackA,
                                        const G4Track &trackB,
                                        const std::vector<G4Track *> *pProducts)
{
  if ((GetMolecule((&trackA))->GetDefinition() == G4Deoxyribose::Definition()) ||
      (GetMolecule((&trackB))->GetDefinition() == G4Deoxyribose::Definition()) ||
      (GetMolecule((&trackA))->GetDefinition() == G4Guanine::Definition()) ||
      (GetMolecule((&trackB))->GetDefinition() == G4Guanine::Definition()) ||
      (GetMolecule((&trackA))->GetDefinition() == G4Thymine::Definition()) ||
      (GetMolecule((&trackB))->GetDefinition() == G4Thymine::Definition()) ||
      (GetMolecule((&trackA))->GetDefinition() == G4Adenine::Definition()) ||
      (GetMolecule((&trackB))->GetDefinition() == G4Adenine::Definition()) ||
      (GetMolecule((&trackA))->GetDefinition() == G4Cytosine::Definition()) ||
      (GetMolecule((&trackB))->GetDefinition() == G4Cytosine::Definition()))
  {

    const G4Track *DNAElement = nullptr;
    const G4Track *radical = nullptr;
    if ((GetMolecule(&trackA)->GetDefinition() == G4Deoxyribose::Definition()) ||
        (GetMolecule(&trackA)->GetDefinition() == G4Guanine::Definition()) ||
        (GetMolecule(&trackA)->GetDefinition() == G4Thymine::Definition()) ||
        (GetMolecule(&trackA)->GetDefinition() == G4Adenine::Definition()) ||
        (GetMolecule(&trackA)->GetDefinition() == G4Cytosine::Definition()))
    {
      DNAElement = &trackA;
      radical = &trackB;
    }
    else
    {
      DNAElement = &trackB;
      radical = &trackA;
    }

    // if (GetMolecule(radical)->GetDefinition() != G4OH::Definition())
    // {
    //   return;
    // }

    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;

    G4ThreeVector localPosDNA = DNAElement->GetPosition();
    G4ThreeVector localPosOH = radical->GetPosition();

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    auto result0 = fPositions0->Nearest(localPosDNA);
    auto result1 = fPositions1->Nearest(localPosDNA);
    auto resultBase0 = fPositionsBase0->Nearest(localPosDNA);
    auto resultBase1 = fPositionsBase1->Nearest(localPosDNA);

    auto resultOH0 = fPositions0->Nearest(localPosOH);
    auto resultOH1 = fPositions1->Nearest(localPosOH);
    auto resultBaseOH0 = fPositionsBase0->Nearest(localPosOH);
    auto resultBaseOH1 = fPositionsBase1->Nearest(localPosOH);
    G4ThreeVector localPos;

    if (((result0)&&(sqrt(result0->GetDistanceSqr()) < 1e-14)) || ((result1)&&(sqrt(result1->GetDistanceSqr()) < 1e-14)) || ((resultBase0)&&(sqrt(resultBase0->GetDistanceSqr()) < 1e-14)) || ((resultBase1)&&(sqrt(resultBase1->GetDistanceSqr()) < 1e-14)))
    {
      localPos = localPosDNA;
    }
    else if ((sqrt(resultOH0->GetDistanceSqr()) < 1e-14) || (sqrt(resultOH1->GetDistanceSqr()) < 1e-14) || (sqrt(resultBaseOH0->GetDistanceSqr()) < 1e-14) || (sqrt(resultBaseOH1->GetDistanceSqr()) < 1e-14))
    {
      localPos = localPosOH;
    }
    else
    {
      G4cout << "position not found " << G4endl;
    }

    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4int step1_eventID = generatorAction->step1_eventID;
    G4int step1_copyNo = generatorAction->step1_copyNo;
    G4int step1_PID = generatorAction->step1_PID;
    G4double step1_time = generatorAction->step1_time;
    G4int step1_primaryID = generatorAction->step1_primaryID;

    analysisManager->FillNtupleIColumn(2, 0, eventID);
    analysisManager->FillNtupleDColumn(2, 1, localPos.x() / nanometer);
    analysisManager->FillNtupleDColumn(2, 2, localPos.y() / nanometer);
    analysisManager->FillNtupleDColumn(2, 3, localPos.z() / nanometer);
    analysisManager->FillNtupleSColumn(2, 4, GetMolecule(DNAElement)->GetName());
    analysisManager->FillNtupleSColumn(2, 5, GetMolecule(radical)->GetName());
    analysisManager->FillNtupleIColumn(2, 6, step1_eventID);
    analysisManager->FillNtupleIColumn(2, 7, step1_copyNo);
    analysisManager->FillNtupleIColumn(2, 8, step1_PID);
    analysisManager->FillNtupleIColumn(2, 9, step1_primaryID);
    analysisManager->FillNtupleDColumn(2, 10, step1_time + DNAElement->GetGlobalTime());
    analysisManager->AddNtupleRow(2);
  }
}

G4Navigator *TimeStepAction::GetNavigator()
{
  static G4ThreadLocal G4Navigator *theNavigator = 0;
  if (!theNavigator)
    theNavigator = new G4Navigator;

  // Make sure current world volume is the one in use
  G4VPhysicalVolume *theWorld =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();

  if (theNavigator->GetWorldVolume() != theWorld)
    theNavigator->SetWorldVolume(theWorld);

  return theNavigator;
}

G4VTouchable *TimeStepAction::CreateTouchableAtPoint(const G4ThreeVector &pos)
{
  G4VTouchable *touchable = new G4TouchableHistory;
  GetNavigator()->LocateGlobalPointAndUpdateTouchable(pos, touchable, false);
  return touchable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
