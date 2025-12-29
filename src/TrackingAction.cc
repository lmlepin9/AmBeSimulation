#include "TrackingAction.hh"

#include <G4StepStatus.hh>
#include <G4IonTable.hh>
#include <G4TrackingManager.hh>

TrackingAction::TrackingAction(EventAction* EventAct, bool NeutronTracking, bool FissFragments, bool ScoreGamma):
  fEventAction(EventAct),
  fNeutronTracking(NeutronTracking),
  fFissFragments(FissFragments),
  fScoreGamma(ScoreGamma)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)//when track is created
{
  if (aTrack->GetParentID()>0)//primaries have parentID==0
  {
    if (fFissFragments)
    {
      if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
        fEventAction->AddFissNeut(aTrack->GetKineticEnergy()/MeV);
      if(aTrack->GetParticleDefinition()->IsGeneralIon() and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
        fEventAction->AddFissIon(aTrack->GetParticleDefinition()->GetAtomicMass());
    }
    if (fNeutronTracking)
    {
      if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()!="nFissionHP")// and aTrack->GetCreatorProcess()->GetProcessName()!="neutronInelastic")
      {
        //G4cout << "TrackingAction: neutron created in volume "<< aTrack->GetVolume()->GetName()<<G4endl;
        //G4cout << "TrackingAction: neutron created in volume "<< aTrack->GetVolume()->GetName()<<" by process "<<aTrack->GetCreatorProcess()->GetProcessName() <<G4endl;
        fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
      }
    }
    if(aTrack->GetParticleDefinition()->GetParticleName()!="neutron" and aTrack->GetParticleDefinition()->GetParticleName()!="gamma")
    {
      fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
    }
    if (fScoreGamma)
    {
      if (aTrack->GetParentID()==2)//Carbon from PrimaryGenerator. Any excited state
      {
        if(aTrack->GetParticleDefinition()->GetParticleName()=="gamma")
        {
          fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
          if (fDBG) G4cout << "G4Tracking:: ParentID->gamma process: "<< aTrack->GetCreatorProcess()->GetProcessName()<<" KE: "<<aTrack->GetKineticEnergy()/keV<<G4endl;
        }
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)//when track is killed
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
