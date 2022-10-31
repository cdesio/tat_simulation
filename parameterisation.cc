

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "parameterisation.hh"

CylinderParameterisation::CylinderParameterisation(G4double total_length,
                                                   G4double inner_radius,
                                                   G4double outer_radius,
                                                   G4int n_div_R,
                                                   G4int n_div_Z,
                                                   G4int n_div_Theta)
{
    this->n_div_R = n_div_R;
    this->n_div_Theta = n_div_Theta;
    this->n_div_Z = n_div_Z;
    NoElements = n_div_R * n_div_Theta * n_div_Z;

    this->map = new G4int *[NoElements];
    G4int row_n = 0;
    for (G4int i = 0; i < n_div_Z; i++)
    {
        for (G4int j = 0; j < n_div_Theta; j++)
        {
            for (G4int k = 0; k < n_div_R; k++)
            {
                this->map[row_n] = new G4int[3];
                this->map[row_n][0] = i;
                this->map[row_n][1] = j;
                this->map[row_n][2] = k;

                // G4cout << map[0] << " " << map[1] << " " << map[2] << G4endl;
                // G4cout << row_n << " " << map[row_n][0] << " " << map[row_n][1] << " " << map[row_n][2] << G4endl;
                row_n += 1;
            }
        }
    }

    // this->vect_trans = new G4int[NoElements];
    // this->vect_dim = new G4int[NoElements];
    // for (G4int i = 0; i < NoElements; i++)
    // {
    //     this->vect_trans[i] = 0;
    //     this->vect_dim[i] = 0;
    // }

    // for (G4int i = 0; i < NoElements; i++)
    // {
    //     for (G4int j = 0; j < 3; j++)
    //     {
    //         G4cout << this->map[i][j] << G4endl;
    //     }
    // }
    // until here.

    this->angle = 360. / n_div_Theta;
    this->half_length = total_length / 2.;
    this->l_n = total_length / n_div_Z;
    this->r_min = inner_radius;
    this->r_max = outer_radius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// CylinderParameterisation::~CylinderParameterisation()
// {
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CylinderParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
    // Note: copyNo will start with zero!
    // G4cout << "computTransformation " << this->vect_trans[copyNo] << G4endl;
    // if (this->vect_trans[copyNo] > this->n_div_Z * this->n_div_Theta)
    //     return;
    // this->vect_trans[copyNo] += 1;
    int i = this->map[copyNo][0];
    int j = this->map[copyNo][1];
    // G4cout << copyNo << " i: " << i << " j: " << j << " " << G4endl;
    G4double Zposition = ((-this->half_length) + ((i + 1) * this->l_n) - this->l_n / 2.);
    // G4cout << copyNo << ", Z: " << Zposition << G4endl;
    G4ThreeVector origin(G4ThreeVector(0., 0., Zposition));
    // G4Rotate3D rotZ(j * angle * deg, G4ThreeVector(0, 0, 1));
    G4RotationMatrix *rotZ = new G4RotationMatrix(); // 0., 0., j * angle * deg);
    rotZ->rotateZ(j * this->angle * deg);
    physVol->SetTranslation(origin);
    physVol->SetRotation(rotZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CylinderParameterisation::ComputeDimensions(G4Tubs &cylinder_element, const G4int copyNo, const G4VPhysicalVolume *) const
{

    // cylinder_element = new G4Tubs("solidTatDetector", 10 * um + (k * 100 / n_x) * um, 10 * um + (100 / n_x) * (k + 1) * um, l_n / 2, 0 * deg, angle * deg);
    //  Note: copyNo will start with zero!
    // G4cout << "computDimensions " << vect_dim[copyNo] << G4endl;
    // if (this->vect_dim[copyNo] > this->n_div_R)
    //     return;
    // this->vect_dim[copyNo] += 1;

    G4int k = this->map[copyNo][2];
    // G4cout << copyNo << ", k: " << k << G4endl;
    G4double inner_radius = this->r_min + (k * ((this->r_max - this->r_min) / this->n_div_R));
    G4double outer_radius = this->r_min + ((this->r_max - this->r_min) / this->n_div_R) * (k + 1);
    // G4cout << inner_radius << " " << outer_radius << G4endl;
    // G4cout << this->r_min << " " << this->r_max << G4endl;

    cylinder_element.SetInnerRadius(inner_radius);
    cylinder_element.SetOuterRadius(outer_radius);
    cylinder_element.SetZHalfLength(this->l_n / 2.);
    cylinder_element.SetStartPhiAngle(0. * deg);
    cylinder_element.SetDeltaPhiAngle(this->angle * deg);
}
// G4double GetInnerRadius();
// G4double GetOuterRadius();
// G4double GetDeltaPhiAngle();