#include "EventAction.hh"
#include "Analysis.hh"

#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4Neutron.hh>

#include <G4SDManager.hh>
#include <G4THitsMap.hh>
#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4Event.hh>

EventAction::EventAction(RunAction* runAction, bool FissFragments, bool NeutronTracking):
  fRunAction(runAction),
  fFissFragments(FissFragments),
  fNeutronTracking(NeutronTracking)
{
  fNeutronDetActivation = false;
  fDBG = false;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::BeginOfEventAction(const G4Event* event)
{
  ClearVars();

  if (event->GetEventID()%10000==0) G4cout << "Event: " << event->GetEventID() << G4endl;
  //G4cout << "Primaries: " << event->GetNumberOfPrimaryVertex() << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  fRunAction->RecordSecondaries(fListSecondaryNeutron,fListSecondaryElectron,fListSecondaryGamma);
  fRunAction->RecordSecondariesEmerging(fListSecondaryNeutronEmerging,fListSecondaryElectronEmerging,fListSecondaryGammaEmerging);

  fRunAction->AddNeutronEmissionSpectrum(fNeutronEmissionSpectrum, fNeutronEmissionSpectrumVer);


  if (fFissFragments)
  {
    FindFission(event->GetEventID()); //find fission fragments - saves to RunAction
    fRunAction->AddFissIon(fFissIon);
    fRunAction->AddFissNeut(fFissNeut);
    fRunAction->AddFissNeutEmerging(fFissNeutEmerging);
  }

  //score emerging neutrons
  if (fEmergingNeutrons.size()!=0) fRunAction->RecordEmergingNeutrons(fEmergingNeutrons);

  fRunAction->FillValidity(fNeutronDetActivation);
  fRunAction->TreeFill(event->GetEventID());
  fNeutronDetActivation = false;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddSecondary(const G4ParticleDefinition* particle, G4double energy)
{
  if (particle == G4Gamma::Definition())
  {
    //fNGammas++;
    //fAverageGammaEnergy += energy;
    fListSecondaryGamma.push_back(energy/MeV);
  }
  else if (particle == G4Electron::Definition())
  {
    //fNElectrons++;
    //fAverageElectronEnergy += energy;
    fListSecondaryElectron.push_back(energy/MeV);
  }
  else if (particle == G4Neutron::Definition())
  {
    //fNNeutrons++;
    //fAverageNeutronEnergy += energy;
    fListSecondaryNeutron.push_back(energy/MeV);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddSecondaryEmerging(const G4ParticleDefinition* particle, G4double energy)
{
  if (particle == G4Gamma::Definition())
  {
    //fNGammas++;
    //fAverageGammaEnergy += energy;
    fListSecondaryGammaEmerging.push_back(energy/MeV);
  }
  else if (particle == G4Electron::Definition())
  {
    //fNElectrons++;
    //fAverageElectronEnergy += energy;
    fListSecondaryElectronEmerging.push_back(energy/MeV);
  }
  else if (particle == G4Neutron::Definition())
  {
    //fNNeutrons++;
    //fAverageNeutronEnergy += energy;
    fListSecondaryNeutronEmerging.push_back(energy/MeV);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ClearVars()
{
  ////Clear vectors
  //secondaries
  fListSecondaryNeutron.clear();
  fListSecondaryNeutronEmerging.clear();
  fListSecondaryElectron.clear();
  fListSecondaryElectronEmerging.clear();
  fListSecondaryGamma.clear();
  fListSecondaryGammaEmerging.clear();
  //Fission process reconstruction
  fTrackList.clear();
  fTrackCreator.clear();
  //Neutron Energy spectrum
  fEmergingNeutrons.clear();

  //Fission
  fFissIon.clear();
  fFissNeut.clear();
  fFissNeutEmerging.clear();

  fNeutronEmissionSpectrum.clear();
  fNeutronEmissionSpectrumVer.clear();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::FindFission(int eventnumber)
{
  int fissionIndex=-1;
  //if (fTrackCreator.size() != fTrackList.size()) return;
  std::vector<int> fissionList;
  //Check particles are fission fragments
  for (unsigned int i=0; i<fTrackCreator.size(); i++)
  {
    if (fTrackCreator[i][0] == G4String("nFissionHP"))
    {
      if (fissionIndex==-1) fissionIndex = i;//first loop
      else if (fissionIndex>=0)
      {
        //G4cout << "DUPLICATE FISSION IN EVENT"<<G4endl;
        fissionList.push_back(fissionIndex);//second loop
      }
    }
  }
  G4String daughter = "";
  int ParentTrackID = -1;
  //first loop
  if (fissionIndex>=0)
  {
    if (fDBG) G4cout << "Event number: " << eventnumber <<G4endl;
    ParentTrackID = fTrackList[fissionIndex][1];
    for (unsigned int i=0; i<fTrackList.size(); i++)
    {
      if (fTrackList[i][1] == ParentTrackID)//particle is fragment
      {
        if (G4String(fTrackCreator[i][1])!=daughter)//skip duplicates
        {
          daughter = G4String(fTrackCreator[i][1]);
          if (fDBG) G4cout <<fTrackCreator[i][0] <<" daughter "<<daughter<<" parent id "<<fissionIndex<<G4endl;
          fRunAction->AddFissionProduct(daughter, fTrackList[i][2]);
        }
      }
    }
  }
  //second loop
  for (unsigned int j=0; j<fissionList.size(); j++)
  {
    ParentTrackID = fTrackList[fissionList[j]][1];
    for (unsigned int i=0; i<fTrackList.size(); i++)
    {
      if (fTrackList[i][1] == ParentTrackID)//particle is fragment
      {
        if (G4String(fTrackCreator[i][1])!=daughter)//skip duplicates
        {
          if (fDBG) G4cout <<fTrackCreator[i][0] <<" daughter "<<fTrackCreator[i][1]<<" parent id "<<fissionList[j]<<G4endl;
          fRunAction->AddFissionProduct(daughter,fTrackList[i][2]);
        }
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::FindParent(int eventnumber)
{
  int fissionIndex=-1;
  //if (fTrackCreator.size() != fTrackList.size()) return;
  std::vector<int> fissionList;
  //Check particles are fission fragments
  for (unsigned int i=0; i<fTrackCreator.size(); i++)
  {
    if (fDBG) G4cout <<fTrackCreator[i][0]<<" "<<fTrackCreator[i][1]<<"\t";
    if (fTrackCreator[i][0] == G4String("photonNuclear") or fTrackCreator[i][0] == G4String("phot") or fTrackCreator[i][0] == G4String("neutronInelastic") or fTrackCreator[i][0] == G4String("eBrem") or fTrackCreator[i][0] == G4String("Radioactivation") or fTrackCreator[i][0] == G4String("hadElastic"))
    {
      if (fissionIndex==-1) fissionIndex = i;//first loop
      else if (fissionIndex>=0)
      {
        G4cout << "DUPLICATE FISSION IN EVENT"<<G4endl;
        fissionList.push_back(fissionIndex);//second loop
      }
    }
  }
  if (fDBG) G4cout<<G4endl;
  G4String daughter = "";
  int ParentTrackID = -1;
  //first loop
  if (fissionIndex>=0)
  {
    if (fDBG) G4cout << "Event number: " << eventnumber <<G4endl;
    ParentTrackID = fTrackList[fissionIndex][1];
    for (unsigned int i=0; i<fTrackList.size(); i++)
    {
      if (fTrackList[i][1] == ParentTrackID)//particle is fragment
      {
        if (G4String(fTrackCreator[i][1])!=daughter)//skip duplicates
        {
          daughter = G4String(fTrackCreator[i][1]);
          if (fDBG) G4cout <<fTrackCreator[i][0] <<" daughter "<<daughter<<" parent id "<<fissionIndex<<G4endl;
          //fRunAction->AddFissionProduct(daughter, fTrackList[i][2]);
        }
      }
    }
  }
  daughter="";
  //second loop
  for (unsigned int j=0; j<fissionList.size(); j++)
  {
    ParentTrackID = fTrackList[fissionList[j]][1];
    for (unsigned int i=0; i<fTrackList.size(); i++)
    {
      if (fTrackList[i][1] == ParentTrackID)//particle is fragment
      {
        if (G4String(fTrackCreator[i][1])!=daughter)//skip duplicates
        {
          if (fDBG) G4cout <<fTrackCreator[i][0] <<" daughter "<<fTrackCreator[i][1]<<" parent id "<<fissionList[j]<<G4endl;
          //fRunAction->AddFissionProduct(daughter,fTrackList[i][2]);
        }
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
