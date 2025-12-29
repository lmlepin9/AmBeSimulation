#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>

class G4LogicalVolume;
class G4GenericMessenger;
class G4Material;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //DetectorConstruction(G4int noWaterBath, G4String IsotopeString);
  DetectorConstruction();

  virtual ~DetectorConstruction();

  void ConstructSDandField() override;

  virtual G4VPhysicalVolume* Construct() override;

  void DetectorMessenger();

  G4Tubs* GetAmBeSolid() {return AbsorberSolid;};
  G4int GetWaterStatus() {return fnoWaterBath;};
  G4VPhysicalVolume* GetWaterTank() {return fWaterTank;};

  void SetWaterBath(G4int noWaterBath) {fnoWaterBath=noWaterBath;};
  void SetIsotope(G4String IsotopeString) {fRadioIsotope=IsotopeString;};
  void SetCasing(G4int CasingSelection) {fCasingSelection=CasingSelection;};
  void SetAzimuthalScoring(G4bool AzimuthalScoring) {fAzimuthalScoring=AzimuthalScoring;};


protected:
    //G4LogicalVolume*  fScoringVolume = nullptr;

private:
  G4GenericMessenger* fMessenger = nullptr;
  G4int fnoWaterBath = 0;
  G4String fRadioIsotope = "241Am";
  G4bool fAzimuthalScoring = false;


  G4bool fDebug;

  G4Tubs* ContainSolid = nullptr;
  G4LogicalVolume* ContainLog = nullptr;
  G4Tubs* AbsorberSolid = nullptr;
  G4LogicalVolume* AbsorberLog = nullptr;

  G4VPhysicalVolume* fWaterTank = nullptr;
  G4int fCasingSelection;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
