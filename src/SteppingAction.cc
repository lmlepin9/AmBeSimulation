#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include <string.h>
#include <G4Step.hh>

#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>
#include <G4LogicalVolume.hh>

#include <G4Event.hh>
#include <G4RunManager.hh>

#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4IonTable.hh>

SteppingAction::SteppingAction(EventAction* eventAction, bool FissFragments, bool NeutronTracking, bool ScoreGamma, bool AzimuthalScoring)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fFissFragments(FissFragments),
  fNeutronTracking(NeutronTracking),
  fScoreGamma(ScoreGamma),
  fAzimuthalScoring(AzimuthalScoring)
{
  fDBG = false;
  fDecayLimit = 30*year;

  fnoWater = fDetectorConstruction->GetWaterStatus();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if (aStep->GetTrack()->GetGlobalTime() > fDecayLimit)//kill step is time above fTimeLimit (30*years)
  {
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);//Disables scoring from any neutron
    return;
  };
  const G4StepPoint* postPoint = aStep->GetPostStepPoint();//EndPoint of step
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();//EndPoint of step
  const G4Track* aTrack = aStep->GetTrack();
  G4LogicalVolume* startVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4String startVolumeName = startVolume->GetName();
  G4double eDepTotal = aStep->GetTotalEnergyDeposit();
  //G4double eDep = aStep->GetEnergyDeposit();
  G4double eDep = eDepTotal;

  const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  if (!aTrack->GetNextVolume())
  {
    //G4cout << "particle out of World"<<G4endl;//DEBUG
    return; //next volume does not exist, i.e. OutOfWorld
  }

  G4LogicalVolume* endVolume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4String endVolumeName = endVolume->GetName();

  G4String particleName = aTrack->GetParticleDefinition()->GetParticleName();

  if (particleName=="neutron")
    if (startVolumeName=="ContainerLogical")
      if (endVolumeName=="WorldLogical")
      {
        fEventAction->ScoreEmergingNeutron(aStep->GetPostStepPoint()->GetKineticEnergy());
      }

  //Tracking emerging secondaries if not Fission
  if (fNeutronTracking)
  {
        if (aTrack->GetParentID()>0)
    {
      if (startVolumeName=="ContainerLogical")
      {
        //if (endVolumeName!="ContainerLogical" or endVolumeName!="AmBeLogical")//WARNING: Don't use, it gives more particles than actually beign created
        if ((fWaterTankPresent and endVolumeName=="waterTankLogical") or (!fWaterTankPresent and endVolumeName=="WorldLogical"))
        {
          if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()!="nFissionHP")
          {
            fEventAction->AddSecondaryEmerging(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
          }
        }
      }
    }
  }
  if (fScoreGamma)
  {
    if (aTrack->GetParentID()>0)//Carbon from PrimaryGenerator. Any excited state
    {
      if (startVolumeName=="ContainerLogical")
      {
        if ((fWaterTankPresent and endVolumeName=="waterTankLogical") or (!fWaterTankPresent and endVolumeName=="WorldLogical"))
        {
          if(aTrack->GetParticleDefinition()->GetParticleName()=="gamma")
          {
            fEventAction->AddSecondaryEmerging(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
          }
        }
      }
    }
  }

  //Fission products check
  if (fFissFragments)
  {
    if (aTrack->GetParentID()>0)
      if (startVolumeName=="ContainerLogical")
        //if (endVolumeName!="ContainerLogical" or endVolumeName!="AmBeLogical")//WARNING: Don't use, it gives more particles than actually beign created
        if (endVolumeName=="WorldLogical")
          if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
            fEventAction->AddFissNeutEmerging(aTrack->GetKineticEnergy()/MeV);

    if (postPoint->GetProcessDefinedStep()->GetProcessName()=="nFissionHP")
    {
      //G4cout << "Fission particle : "<< aTrack->GetParticleDefinition()->GetParticleName() <<G4endl;
      if (fDBG) G4cout << aTrack->GetParticleDefinition()->GetParticleName()<<"Process name: "<<postPoint->GetProcessDefinedStep()->GetProcessName()<<G4endl;
    }
    if(aTrack->GetParticleDefinition()->IsGeneralIon())
    {
      if ((aTrack->GetParticleDefinition()->GetAtomicMass()==12 and aTrack->GetParticleDefinition()->GetAtomicNumber()==6) or (aTrack->GetParticleDefinition()->GetAtomicMass()==16 and aTrack->GetParticleDefinition()->GetAtomicNumber()==8))
      {}
      else
      {
        //G4cout << aTrack->GetParticleDefinition()->GetParticleName()<<"Process name: "<<postPoint->GetProcessDefinedStep()->GetProcessName()<<G4endl;
        G4int parentID = aTrack->GetParentID();
        G4int trackID = aTrack->GetTrackID();
        fEventAction->AddTrack(parentID,trackID,aTrack->GetCreatorProcess()->GetProcessName(),aTrack->GetParticleDefinition()->GetAtomicMass(),aTrack->GetParticleDefinition()->GetParticleName());
      }
    }
  }

  //Neutron Spectrum scoring Perpendicular and Vertical
  if (fAzimuthalScoring)
  {
    if (particleName=="neutron")
    {
      if (((fWaterTankPresent and startVolumeName=="waterTankLogical") or (!fWaterTankPresent and startVolumeName=="WorldLogical")) and endVolumeName=="EnerSphereLogical")
      {
        //G4cout << "SteppingAction: Azimuthal Spectrum: scoring neutron"<<G4endl;
        fEventAction->AddNeutronEmissionSpectrum(aTrack->GetKineticEnergy(), aTrack->GetPosition());
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
