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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: Henri Payno and Yann Perrot
//
//
/// \file ClusterSBPoints.cc
/// \brief Implementation of the ClusterSBPoints class

#include "ClusterSBPoints_.hh"
#include <iostream>
#include <numeric>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusterSBPoints::ClusterSBPoints(std::set<SBPoint *> pSBPoints, bool pContinuous) : fpRegisteredSBPoints(), fContinuous(pContinuous)
{
  std::set<SBPoint *>::iterator itPt;
  for (itPt = pSBPoints.begin(); itPt != pSBPoints.end(); ++itPt)
  {
    AddSBPoint(*itPt);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusterSBPoints::~ClusterSBPoints()
{
  Clear();
}

void ClusterSBPoints::Clear()
{
  std::set<SBPoint *>::iterator itPt;
  for (itPt = fpRegisteredSBPoints.begin();
       itPt != fpRegisteredSBPoints.end();
       ++itPt)
  {
    (*itPt)->CleanCluster();
  }
  fpRegisteredSBPoints.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusterSBPoints::AddSBPoint(SBPoint *pSBPoint)
{
  assert(pSBPoint);
  fpRegisteredSBPoints.insert(pSBPoint);
  pSBPoint->SetCluster(this);

  UpdateDoubleStrand();
  UpdateStrandSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t ClusterSBPoints::GetBarycenter() const
{
  int64_t copyNum = 0;

  std::set<SBPoint *>::iterator itSDSPt;
  for (itSDSPt = fpRegisteredSBPoints.begin();
       itSDSPt != fpRegisteredSBPoints.end();
       ++itSDSPt)
  {
    copyNum += (*itSDSPt)->GetPosition();
  }

  return copyNum / fpRegisteredSBPoints.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ClusterSBPoints::UpdateStrandSource()
{
  std::set<SBPoint *>::iterator itSDSPt;

  std::vector<int> sourcesStrand0;
  std::vector<int> sourcesStrand1;

  int i = 0;
  for (itSDSPt = fpRegisteredSBPoints.begin();
       itSDSPt != fpRegisteredSBPoints.end();
       ++itSDSPt)
  {
    if ((*itSDSPt)->GetStrandSource() == 0)
      std::cout << "ERROR" << std::endl;

    if ((*itSDSPt)->GetTouchedStrand() == 0)
    {
      sourcesStrand0.push_back((*itSDSPt)->GetStrandSource());
    }
    else
    {
      sourcesStrand1.push_back((*itSDSPt)->GetStrandSource());
    }
    i++;
  }

  double meanSource0 = 1.0 * std::accumulate(sourcesStrand0.begin(), sourcesStrand0.end(), 0) / sourcesStrand0.size();
  double meanSource1 = 1.0 * std::accumulate(sourcesStrand1.begin(), sourcesStrand1.end(), 0) / sourcesStrand1.size();

  if ((sourcesStrand0.size() == 0) || (sourcesStrand1.size() == 0)) // is a complex SSB
  {
    if ((meanSource0 == 1.0) || ((meanSource1 == 1.0))) // direct
      fSource = 1;
    else if ((meanSource0 == 2.0) || ((meanSource1 == 2.0))) // indirect
      fSource = 2;
    else
      fSource = 4;
  }
  else // is DSB
  {
    if ((meanSource0 == 1.0) && ((meanSource1 == 1.0))) // direct
      fSource = 1;
    else if ((meanSource0 == 2) && ((meanSource1 == 2))) // indirect
      fSource = 2;
    else if (((meanSource0 == 2) && ((meanSource1 == 1)) || ((meanSource0 == 1) && ((meanSource1 == 2))))) // hybrid
      fSource = 3;
    else
      fSource = 4;
  }

  return;
}

void ClusterSBPoints::UpdateDoubleStrand()
{
  fIsDoubleSB = false;
  bool firstStrandTouch = false;
  bool secondStrandTouch = false;

  std::set<SBPoint *>::iterator itSDSPt;
  for (itSDSPt = fpRegisteredSBPoints.begin();
       itSDSPt != fpRegisteredSBPoints.end();
       ++itSDSPt)
  {
    // if the SDSPoint is localized on the first strand
    if (((*itSDSPt)->GetTouchedStrand() == 0) && !firstStrandTouch)
    {
      firstStrandTouch = true;
      if (secondStrandTouch)
      {
        fIsDoubleSB = true;
        return;
      }
    }
    // if the SDSPoint is localized on the second strand
    if (((*itSDSPt)->GetTouchedStrand() == 1) && !secondStrandTouch)
    {
      secondStrandTouch = true;
      if (firstStrandTouch)
      {
        fIsDoubleSB = true;
        return;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool AreOnTheSameCluster(const SBPoint *pPt1, const SBPoint *pPt2,
                         int64_t pMinBP, bool fContinuous)
{
  assert(pPt1);
  assert(pPt2);

  int64_t copy1 = pPt1->GetPosition();
  int64_t copy2 = pPt2->GetPosition();

  if (fContinuous == false)
  {
    // check the BP are on the same segments of chromatin fibre, clustered damage does not wrap between segments
    int64_t chromatinSegment1 = copy1 / 21168;
    int64_t chromatinSegment2 = copy2 / 21168;

    // check the BP are on the same nucleosome, clustered damage does not wrap between nucleosomes due to the large gap because linking bp are not included
    int64_t nucleosome1 = copy1 / 147;
    int64_t nucleosome2 = copy2 / 147;

    if ((abs(copy1 - copy2) <= pMinBP) && (chromatinSegment1 == chromatinSegment2) && (nucleosome1 == nucleosome2))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    if (abs(copy1 - copy2) <= pMinBP)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusterSBPoints::FindAllPointsPossible(std::vector<SBPoint *> *pPtsToCheck,
                                            int64_t pMinPts, double pMinDist)
{
  assert((uint64_t)pMinPts > this->GetSize());
  std::vector<SBPoint *>::iterator itPt = pPtsToCheck->begin();
  while (itPt != pPtsToCheck->end())
  {

    // If 1- each SBpoint is part of only one cluster
    //    2- the point isn't already in the cluster
    //    3- the point is close enough of the barycenter
    if ((!(*itPt)->HasCluster()) && (fpRegisteredSBPoints.find(*itPt) == fpRegisteredSBPoints.end()) && HasInBarycenter(*itPt, pMinDist)) // first version used HasIn method
    {
      // the point is added
      this->AddSBPoint(*itPt);
      if (this->GetSize() >= (uint64_t)pMinPts)
      {
        return;
      }
      // restart from scratch
      itPt = pPtsToCheck->begin();
    }
    else
    {
      ++itPt;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusterSBPoints::HasIn(const SBPoint *pPtToCheck, double pMinDist)
{
  // check if the given point is near one of the cluster's point
  std::set<SBPoint *>::iterator itClusPt;
  for (itClusPt = fpRegisteredSBPoints.begin();
       itClusPt != fpRegisteredSBPoints.end();
       ++itClusPt)
  {
    // if are two different pts
    if ((*pPtToCheck != *(*itClusPt)))
    {
      // if close enought of an include point of the cluster
      if (AreOnTheSameCluster(pPtToCheck, *itClusPt, pMinDist, fContinuous))
      {
        return true;
      }
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusterSBPoints::HasInBarycenter(const SBPoint *pPtToCheck,
                                      double pMinBP)
{

  int64_t copy1 = pPtToCheck->GetPosition();
  int64_t copy2 = this->GetBarycenter();

  if (abs(copy1 - copy2) <= pMinBP)

  {
    return true;
  }
  else
  {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// this will insert all registredSBPoint
/// from the given cluster to this cluster.
void ClusterSBPoints::MergeWith(ClusterSBPoints *pCluster)
{
  std::set<SBPoint *> points = pCluster->GetRegistredSBPoints();
  pCluster->Clear();
  std::set<SBPoint *>::iterator itPt;
  for (itPt = points.begin(); itPt != points.end(); ++itPt)
  {
    this->AddSBPoint(*itPt);
  }
}