#ifndef PARAMETERISATION_HH
#define PARAMETERISATION_HH

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

class CylinderParameterisation : public G4VPVParameterisation
{
public:
    CylinderParameterisation(G4double total_length,
                             G4double inner_radius,
                             G4double outer_radius,
                             G4int n_div_R,
                             G4int n_div_Z,
                             G4int n_div_Theta);

    ~CylinderParameterisation()
    {
    }

    void ComputeTransformation(const G4int copyNo,
                               G4VPhysicalVolume *physVol) const override;

    void ComputeDimensions(G4Tubs &trackerLayer, const G4int copyNo,
                           const G4VPhysicalVolume *physVol) const override;

    // G4double GetInnerRadius();
    // G4double GetOuterRadius();
    // G4double GetDeltaPhiAngle();

private:
    G4double l_n;   //  length of z-element
    G4double angle; // k-th angle
    G4double half_length;
    G4int NoElements;
    G4int **map;
    // G4int *vect_trans, *vect_dim;
    G4double r_min;
    G4double r_max;
    int n_div_R;
    int n_div_Z;
    int n_div_Theta;
};

#endif
