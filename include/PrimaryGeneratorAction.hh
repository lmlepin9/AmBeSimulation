#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4GeneralParticleSource.hh>
#include <globals.hh>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <G4SystemOfUnits.hh>
#include <G4GenericMessenger.hh>

#include "RunAction.hh"
#include "NNDCLoader.hh"
#include "DetectorConstruction.hh"

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class G4Box;
class G4GenericMessenger;
class RunAction;
class NNDCLoader;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(RunAction* runAction, DetectorConstruction* detectorConstruction, G4String IsotopeString, G4int NumberOfThreads);
  ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event* anEvent) override;

  const G4ParticleGun *GetParticleGun() const { return fParticleGun; }
  //const G4GeneralParticleSource* GetParticleGun() const { return fGPS; }


  void PrimaryGeneratorMessenger();

private:
  G4GenericMessenger* fMessenger = nullptr;
  G4ParticleGun *fParticleGun  = nullptr;
  G4Box *fEnvelopeBox  = nullptr;
  G4int fThreadid = -1;

  Double_t GetNeutron(Double_t alphaEnergy);
  Double_t InteractionE(Double_t Eint);
  Double_t InteractionTheta(Double_t E_inc);
  Double_t Eloss(Double_t Ea, Double_t z);
  Double_t ConvertCMtoLab(Double_t Eint,Double_t thetaCM);
  Double_t NeutronEnergy(Double_t Eint, Double_t thetaCM);

  /*
  G4GeneralParticleSource *fGPS = nullptr;
  G4SingleParticleSource *fCarbonGPS = nullptr;
  G4SingleParticleSource *fNeutronGPS = nullptr;
  G4SingleParticleSource *fGamma_59_54, *fGamma_26_34, *fGamma_13_9 = nullptr;
  //*/
  std::vector <G4GeneralParticleSource*> fGPS;
  std::vector <G4SingleParticleSource*> fCarbonGPS;
  std::vector <G4SingleParticleSource*> fNeutronGPS;
  std::vector <G4SingleParticleSource*> fGamma_59_54;
  std::vector <G4SingleParticleSource*> fGamma_26_34;
  std::vector <G4SingleParticleSource*> fGamma_13_9;

  G4bool f241AmGammaEnabled = false;
  G4String fRadioIsotope = "241Am";
  G4double singleNeutronEnergy = -1;
  G4int fNumberOfThreads = -1;

  DetectorConstruction *fDetectorConstruction = nullptr;
  RunAction *fRunAction = nullptr;
  NNDCLoader *fNNDCAlpha = nullptr;

  G4int fDBG {0};

  Bool_t secondks=false;
  Bool_t excited=false;
  Bool_t excited2=false;
  Bool_t excited3=false;

  TGraph *gXS_0 = nullptr;
  TGraph *gXS_1 = nullptr;
  TGraph *gXS_2 = nullptr;
  TGraph *gXS_3 = nullptr;
  TGraph *gstopping = nullptr;
  TGraph *gXS_t = nullptr;
  //Tgraph2DErrors *g_gs;
  //TGraph2DErrors *gp;
  //TGraph2DErrors *gpp;
  std::vector<TGraph*> g_gs;
  std::vector<TGraph*> gp;
  std::vector<TGraph*> gpp;
  std::vector<TGraph*> gppp;
  Double_t fIntXSMax {-1};

  //Constants
  Double_t d2r=atan(1.)/45.;
  Double_t pi=4.*atan(1.);
  Double_t amu=931.4941;//MeV/u

  //12C data
  G4int Z = 6;
  G4int A = 12;
  G4double ionCharge = 0.*eplus;
  G4double C_Excit_GND = 0.*MeV;
  G4double C_Excit_1ST = 4.43982*MeV;
  G4double C_Excit_2ND = 7.65407*MeV;
  G4double C_Excit_3RD = 9.641*MeV;
  Double_t C_Excit_GND_MeV = 0.;
  Double_t C_Excit_1ST_MeV = 4.43982;
  Double_t C_Excit_2ND_MeV = 7.65407;
  Double_t C_Excit_3RD_MeV = 9.641;

  Double_t m_4He=4.002603;//u -- Mass 4He
  Double_t m_9Be=9.012182;//u -- Mass 9Be
  Double_t m_n=1.0008665;//u -- Mass neutron
  Double_t m_12C=12.;//u -- Mass 12C
  Double_t Q_9BeAN=+5.702049;//MeV

  Double_t max_depth = 0.03704365;// old value: 0.03639234;//mm according to LISE++ at E=Einitial for 241Am:110=15.999O:220=9.012Be:890 and density 2.9265704*g/cm3
  Double_t Einitial=5.48556;//MeV
  Double_t Eth=0. ;//Threshold energy

  G4double alphaPerNeutron = 35890.0/2.27;
  G4double r_max, z_max;

  TF1 *fLeg8 = nullptr;

  G4double USF()
    {
      Bool_t trials = true;
      Int_t loopcounter=0;
      while(trials) //Randomly choose a depth into the target (rather than energy to keep ame target thickness)
      {
        G4double energy = G4UniformRand()*10;
        G4double eneProb = std::sinh(sqrt(2.249*energy))*std::exp(-energy/0.988);
        //G4double eneProb = std::sinh(sqrt(2*energy))*std::exp(-energy); //constant 4.75*
        //G4double eneProb = (1.088*std::sinh(sqrt(1.92*energy))+2.937*std::sinh(std::sqrt(1.28*energy)))*std::exp(-energy/0.875); //normalisation 1/4.025
        if (eneProb>=G4UniformRand()*1.5)
        {
          return energy*MeV;
          trials = false;
        }
        if (loopcounter>1000) return -1;
        loopcounter++;
      }
      return -1;
    }

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
