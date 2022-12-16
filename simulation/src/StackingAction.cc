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
/// \file StackingAction.cc
/// \brief Implementation of the StackingAction class

#include "StackingAction.hh"
#include "G4DNAChemistryManager.hh"
#include "G4StackManager.hh"
#include "G4ITTransportationManager.hh"
#include "G4ITTrackHolder.hh"
#include "DetectorConstruction.hh"
#include "G4DNAMolecule.hh"
#include "G4Molecule.hh"
#include "G4KDTree.hh"
#include "G4KDTreeResult.hh"
// #include <G4UnitsTable.hh>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction *pDetector)
    : G4UserStackingAction(), fpDetector(pDetector)
{
    fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
    dkill = ((DetectorConstruction *)fpDetector)->Getdkill();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::NewStage()
{

    if (stackManager->GetNTotalTrack() == 0)
    {
        ///  deoxyribose & bases
        auto fPositions0 = ((DetectorConstruction *)fpDetector)->fPositions0;
        auto fPositionsBase0 = ((DetectorConstruction *)fpDetector)->fPositionsBase0;
        auto fPositions1 = ((DetectorConstruction *)fpDetector)->fPositions1;
        auto fPositionsBase1 = ((DetectorConstruction *)fpDetector)->fPositionsBase1;

        G4KDTree *tree0 = fpEventAction->fPositions0Event;
        G4KDTree *tree1 = fpEventAction->fPositions1Event;
        G4KDTree *treeBase0 = fpEventAction->fPositionsBase0Event;
        G4KDTree *treeBase1 = fpEventAction->fPositionsBase1Event;

        std::vector<G4ThreeVector> forTree0;
        std::vector<G4ThreeVector> forTree1;
        std::vector<G4ThreeVector> forTree0Base;
        std::vector<G4ThreeVector> forTree1Base;

        auto delayList = G4ITTrackHolder::Instance()->GetDelayedLists();
        for (auto &delayedmap_it : delayList)
        {
            for (auto &trackList : delayedmap_it.second)
            {
                if (nullptr == trackList.second)
                {
                    return;
                }
                G4TrackList::iterator itt = trackList.second->begin();
                G4TrackList::iterator endd = trackList.second->end();
                for (; itt != endd; ++itt)
                {
                    G4Track *track = *itt;
                    G4ThreeVector pos = track->GetPosition();

                    G4KDTreeResultHandle result0 = fPositions0->NearestInRange(pos, dkill); // find all sugars within dkill of the track
                    G4KDTreeResultHandle result1 = fPositions1->NearestInRange(pos, dkill); // find all sugars within dkill of the track

                    // Remove all radicals more than dkill from DNA
                    if ((result0->GetSize() == 0) && (result1->GetSize() == 0))
                    {
                        track->SetTrackStatus(fKillTrackAndSecondaries);
                    }
                    else
                    {
                        for (result0->Rewind(); !result0->End(); result0->Next())
                        {
                            auto node = result0->GetNode();
                            auto p = G4ThreeVector((*node)[0], (*node)[1], (*node)[2]);

                            if (std::find(forTree0.begin(), forTree0.end(), p) == forTree0.end()) // if sugar not already in vector
                            {
                                forTree0.push_back(p);
                            }
                        }
                        for (result1->Rewind(); !result1->End(); result1->Next())
                        {
                            auto node = result1->GetNode();
                            auto p = G4ThreeVector((*node)[0], (*node)[1], (*node)[2]);

                            if (std::find(forTree1.begin(), forTree1.end(), p) == forTree1.end()) // if sugar not already in vector
                            {
                                forTree1.push_back(p);
                            }
                        }
                        G4KDTreeResultHandle resultBase0 = fPositionsBase0->NearestInRange(pos, dkill); 
                        G4KDTreeResultHandle resultBase1 = fPositionsBase1->NearestInRange(pos, dkill); 
                        for (resultBase0->Rewind(); !resultBase0->End(); resultBase0->Next())
                        {
                            auto node = resultBase0->GetNode();
                            auto p = G4ThreeVector((*node)[0], (*node)[1], (*node)[2]);

                            if (std::find(forTree0Base.begin(), forTree0Base.end(), p) == forTree0Base.end()) // if sugar not already in vector
                            {
                                forTree0Base.push_back(p);
                            }
                        }
                        for (resultBase1->Rewind(); !resultBase1->End(); resultBase1->Next())
                        {
                            auto node = resultBase1->GetNode();
                            auto p = G4ThreeVector((*node)[0], (*node)[1], (*node)[2]);

                            if (std::find(forTree1Base.begin(), forTree1Base.end(), p) == forTree1Base.end()) // if sugar not already in vector
                            {
                                forTree1Base.push_back(p);
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < forTree0.size(); ++i)
        {
                tree0->Insert(forTree0[i]);

            G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Deoxyribose::Definition())),
                                                            1 * CLHEP::picosecond, forTree0[i], -100);
        }
        for (int i = 0; i < forTree1.size(); ++i)
        {
                tree1->Insert(forTree1[i]);

            G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Deoxyribose::Definition())),
                                                            1 * CLHEP::picosecond, forTree1[i], -100);
        }
        for (int i = 0; i < forTree0Base.size(); ++i)
        {
            treeBase0->Insert(forTree0Base[i]);
            G4double R = G4UniformRand();
            if (R < 0.5) // Add bases randomly
            {
                G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Cytosine::Definition())),
                                                                1 * CLHEP::picosecond, forTree0Base[i], -100);
            }
            else
            {
                G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Adenine::Definition())),
                                                                1 * CLHEP::picosecond, forTree0Base[i], -100);
            }
        }
        for (int i = 0; i < forTree1Base.size(); ++i)
        {
            treeBase1->Insert(forTree1Base[i]);
            G4double R = G4UniformRand();
            if (R < 0.5) // Add bases randomly
            {
                G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Guanine::Definition())),
                                                                1 * CLHEP::picosecond, forTree1Base[i], -100);
            }
            else
            {
                G4DNAChemistryManager::Instance()->PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Thymine::Definition())),
                                                                1 * CLHEP::picosecond, forTree1Base[i], -100);
            }
        }

        G4ITTransportationManager::GetTransportationManager()->SetWorldForTracking(
            G4ITTransportationManager::GetTransportationManager()->GetParallelWorld("ChemistryWorld"));
        G4DNAChemistryManager::Instance()->Run(); // starts chemistry
    }
}
