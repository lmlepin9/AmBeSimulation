#ifndef TRACKINGACTION_HH
#define TRACKINGACTION_HH

#include "EventAction.hh"

#include <G4UserTrackingAction.hh>
#include <G4Track.hh>

#include <globals.hh>

class EventAction;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

public:
  TrackingAction(EventAction* EventAct, bool FissFragments, bool ScoreGamma);
  ~TrackingAction() override = default;

  void PreUserTrackingAction(const G4Track* aTrack) override;
  void PostUserTrackingAction(const G4Track* aTrack) override;
  std::map<G4int, G4int> parentInfo; // Maps track ID to parent ID

private:
  EventAction* fEventAction = nullptr;
  G4bool fDBG = false;
  G4bool fFissFragments = false;
  G4bool fScoreGamma = false;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
