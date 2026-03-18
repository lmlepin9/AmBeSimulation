#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include <G4VUserActionInitialization.hh>
#include <globals.hh>
#include "DetectorConstruction.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(int rank, int NumberOfThreads, DetectorConstruction* detectorConstruction, bool FissFragments, bool NeutronTracking, G4String RadioIsotope, bool InitialNeutrons, bool ScoreGamma, bool AzimuthalScoring, bool SaveEmerging);
  ~ActionInitialization();
  //virtual ~ActionInitialization();

  void Build() const override;

  void BuildForMaster() const override;

private:
  G4int fNumberOfThreads { -1 };
  G4bool fFissFragments=false;
  G4bool fNeutronTracking=false;
  G4int fRank = 0;
  G4String fRadioIsotope;
  G4bool fInitialNeutrons;
  G4bool fScoreGamma=false;
  G4bool fAzimuthalScoring=false;
  G4bool fSaveEmerging=false;


  DetectorConstruction *fDetectorConstruction = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
