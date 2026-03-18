#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "RunAction.hh"

#include <globals.hh>
#include <vector>
#include <math.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include <G4UserEventAction.hh>
#include <G4ThreeVector.hh>


#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4Positron.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4Proton.hh>
#include <G4PionMinus.hh>
#include <G4PionPlus.hh>
#include <G4Alpha.hh>
#include <G4GenericIon.hh>
#include <G4He3.hh>
#include <G4Deuteron.hh>
#include <G4Triton.hh>
#include <G4Neutron.hh>



class RunAction;
class Cartesian3D;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction* runAction, bool NeutronTracking, bool FissFragments);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  RunAction* GetRunAction() {return fRunAction;};


  //void SetPrimaryVertexEnergy(G4double PrimaryVertexEnergy);
  //secondaries
  void AddSecondary(const G4ParticleDefinition*, G4double energy);
  void AddSecondaryEmerging(const G4ParticleDefinition*, G4double energy);

  //Fission process reconstruction
  void AddTrack(G4int parent, G4int track, G4String creatorproc, G4int atomicMass, G4String particlename)
    {
      fTrackList.push_back({parent,track,atomicMass});
      fTrackCreator.push_back({creatorproc,particlename});
    };

  void FindFission(int eventnumber);
  void FindParent(int eventnumber);

  void AddFissIon(G4int mass){ fFissIon.push_back(mass); };
  void AddFissNeut(G4double ene){ fFissNeut.push_back(ene); };
  void AddFissNeutEmerging(G4double ene){ fFissNeutEmerging.push_back(ene); };
  void AddNeutronEmissionSpectrum(G4double KE, G4ThreeVector ver)
    {
      fNeutronEmissionSpectrum.push_back((Double_t)(KE/MeV));
      fNeutronEmissionSpectrumVer.push_back((ROOT::Math::XYZVectorD)ver/mm);
    };


  //Neutron Energy spectrum
  void ScoreEmergingNeutron(G4double KE)
    {fEmergingNeutrons.push_back((Double_t)(KE/MeV));};

  //Tree Utils
  void ClearVars();

private:
  RunAction* fRunAction = nullptr;
  //settables
  G4bool fDBG;
  //parameters
  G4bool fFissFragments = false;
  G4bool fNeutronTracking = false;

  G4bool fNeutronDetActivation;


  Double_t gammaEnergy=0.;

  ////Branches
  //secondaries
  std::vector<Double_t> fListSecondaryNeutron;
  std::vector<Double_t> fListSecondaryNeutronEmerging;
  std::vector<Double_t> fListSecondaryElectron;
  std::vector<Double_t> fListSecondaryElectronEmerging;
  std::vector<Double_t> fListSecondaryGamma;
  std::vector<Double_t> fListSecondaryGammaEmerging;
  //Fission process reconstruction
  std::vector<std::array<G4int, 3>> fTrackList;
  std::vector<std::array<G4String, 2>> fTrackCreator;
  //Neutron Energy spectrum
  std::vector<Double_t> fEmergingNeutrons;

  std::vector<Int_t> fFissIon;
  std::vector<Double_t> fFissNeut;
  std::vector<Double_t> fFissNeutEmerging;

  std::vector<Double_t> fNeutronEmissionSpectrum;
  std::vector<ROOT::Math::XYZVectorD> fNeutronEmissionSpectrumVer;

  // Emerging particle record
  std::vector<Int_t> fEmergingId;
  std::vector<Int_t> fEmergingParentId;
  std::vector<Int_t> fEmergingPDG;
  std::vector<ROOT::Math::XYZTVector> fEmergingPos;
  std::vector<ROOT::Math::XYZTVector> fEmergingP; 
  std::vector<std::string> fEmergingProcess; 

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
