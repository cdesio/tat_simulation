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
/// \file SBPoint.hh
/// \brief Definition of the SBPoint class

#ifndef SB_POINT_HH
#define SB_POINT_HH

#include <assert.h>
#include <stdint.h>

#include <cinttypes>

class ClusterSBPoints;
/// \brief defines a point of energy deposition which defines a damage to the DNA.
class SBPoint
{
public:
  /// \brief constructor
  SBPoint(unsigned int, int64_t pPos, int64_t strand, int64_t source );
  /// \brief destructor
  ~SBPoint();

  // Get methods
  int64_t GetID() const
  {
    return fId;
  }
  int64_t GetPosition() const {
    return fPosition;
  }
  ClusterSBPoints* GetCluster() const {
    return fpCluster;
  }
  int64_t GetTouchedStrand() const {
    return fStrand;
  }

  // Set methods
  void SetCluster(ClusterSBPoints* pCluster)
  {
    assert(pCluster); fpCluster = pCluster;
  }

  void SetPosition(int64_t pPos)
  {
    fPosition=pPos;
  }

  bool HasCluster() const {
    return fpCluster != 0;
  }
  void CleanCluster() {
    fpCluster = 0;
  }
  int64_t GetStrandSource(){return fSBsource;}

  bool operator != (const SBPoint& ) const;
  bool operator == (const SBPoint& ) const;
  bool operator < (const SBPoint& ) const;
  bool operator > (const SBPoint& ) const;

private:

  uint64_t fId;             //ID
  int64_t fPosition;      //Position (copy number)
  ClusterSBPoints* fpCluster;   // Associated clustered points
  int64_t fStrand;                // Strand
  int64_t fSBsource{0}; // 0 unknown, 1 direct, 2 indirect


};

#endif // SB_POINT_HH