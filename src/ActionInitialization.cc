#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"

ActionInitialization::ActionInitialization(int rank, int NumberOfThreads, DetectorConstruction* detectorConstruction, bool FissFragments, bool NeutronTracking, G4String RadioIsotope, bool InitialNeutrons, bool ScoreGamma, bool AzimuthalScoring, bool SaveEmerging):
  G4VUserActionInitialization(),
  fRank(rank),
  fNumberOfThreads(NumberOfThreads),
  fDetectorConstruction(detectorConstruction),
  fFissFragments(FissFragments),
  fNeutronTracking(NeutronTracking),
  fRadioIsotope(RadioIsotope),
  fInitialNeutrons(InitialNeutrons),
  fScoreGamma(ScoreGamma),
  fAzimuthalScoring(AzimuthalScoring),
  fSaveEmerging(SaveEmerging)

{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{

  RunAction* theRunAction = new RunAction(fRank, fNumberOfThreads, fFissFragments, fNeutronTracking, fInitialNeutrons, fRadioIsotope, fScoreGamma, fAzimuthalScoring,fSaveEmerging);
  SetUserAction(theRunAction);

  PrimaryGeneratorAction* thePrimaryGenerator = new PrimaryGeneratorAction(theRunAction, fDetectorConstruction, fRadioIsotope, fNumberOfThreads);
  SetUserAction(thePrimaryGenerator);

  EventAction* theEventAction = new EventAction(theRunAction, fFissFragments, fNeutronTracking);
  SetUserAction(theEventAction);

  TrackingAction* theTrackingAction = new TrackingAction(theEventAction, fNeutronTracking, fFissFragments, fScoreGamma);
  SetUserAction(theTrackingAction);

  SteppingAction* theSteppingAction = new SteppingAction(theEventAction, fDetectorConstruction, fFissFragments, fNeutronTracking, fScoreGamma, fAzimuthalScoring);
  SetUserAction(theSteppingAction);

  StackingAction* theStackingAction = new StackingAction(theRunAction,theEventAction);
  SetUserAction(theStackingAction);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  // By default, don't do anything. This applies only in MT mode:
  SetUserAction(new RunAction(fRank, fNumberOfThreads, fFissFragments, fNeutronTracking, fInitialNeutrons, fRadioIsotope, fScoreGamma, fAzimuthalScoring,fSaveEmerging));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
