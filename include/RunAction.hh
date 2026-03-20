#ifndef RUNACTION_HH
#define RUNACTION_HH

#include <G4UserRunAction.hh>
#include <G4Run.hh>
#include <G4ParticleDefinition.hh>
#include <G4Accumulable.hh>
#include <G4GenericMessenger.hh>

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

#include <vector>
#include <globals.hh>

#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TGraph2DErrors.h>

class Cartesian3D;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:
  RunAction(int rank, int fNumberOfThreads, bool FissFragments, bool NeutronTracking, bool InitialNeutrons, G4String RadioIsotope, bool ScoreGamma, bool AzimuthalScoring, bool SaveEmerging);


  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void RunActionMessenger();

  void FillInitialNeutron(Double_t EnerNeutron)
    {fbranchEInitN = EnerNeutron;};
  void RecordSecondaries(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma);
  void RecordSecondariesEmerging(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma);


  void RecordEmerging(std::vector<Int_t> EmergingEventId,
                      std::vector<Int_t> EmergingId,
                      std::vector<Int_t> EmergingParentId,
                      std::vector<Int_t> EmergingPDG,
                      std::vector<ROOT::Math::XYZTVector> EmergingPos,
                      std::vector<ROOT::Math::XYZTVector> EmergingP,
                      std::vector<std::string> EmergingProcess);

  //Neutron Energy spectrum
  void RecordEmergingNeutrons(std::vector<Double_t> emergingNeutrons)
    {fEmergingNeutrons=emergingNeutrons;};
  //Fission Products
  void AddFissionProduct(G4String daughter, G4int mass)
    {
      fbranchFissionProduct.push_back((std::string)daughter);
      fbranchFissionProductMass.push_back(mass);
    };

  void AddFissIon(std::vector<Int_t> mass){ fbranchFissIon = mass; };
  void AddFissNeut(std::vector<Double_t> ene){ fbranchFissNeut = ene; };
  void AddFissNeutEmerging(std::vector<Double_t> ene){ fbranchFissNeutEmerging = ene; };
  void AddNeutronEmissionSpectrum(std::vector<Double_t> E, std::vector<ROOT::Math::XYZVectorD> ver)
    {
      fbranchNEmissionSpec = E;
      fbranchNEmissionSpecVer = ver;
    };


  //TreeUtils
  void TreeFill(G4int eventNumber);
  void ClearBranches();
  void FillValidity(Bool_t status)
    {fbranchValid = status;};


  TGraph* GetgXS_0() {return gXS_0;};
  TGraph* GetgXS_1() {return gXS_1;};
  TGraph* GetgXS_2() {return gXS_2;};
  TGraph* GetgXS_3() {return gXS_3;};
  TGraph* Getgstopping() {return gstopping;};
  TGraph* GetgXS_t() {return gXS_t;};
  //TGraph2DErrors* Getg_gs() {return g_gs;};
  //TGraph2DErrors* Getgp() {return gp;};
  //TGraph2DErrors* Getgpp() {return gpp;};
  std::vector<TGraph*> Getg_gs() {return g_gs;};
  std::vector<TGraph*> Getgp() {return gp;};
  std::vector<TGraph*> Getgpp() {return gpp;};
  std::vector<TGraph*> Getgppp() {return gppp;};

private:
  G4int fDBG = 0;
  G4GenericMessenger* fMessenger = nullptr;

  G4int fNumberOfThreads {-1};
  G4bool fFissFragments = false;
  G4bool fNeutronTracking = false;
  G4int fRank = 0;
  G4bool fInitialNeutrons = false;
  G4String fRadioIsotope;
  G4bool fScoreGamma = false;
  G4bool fAzimuthalScoring = false;
  G4bool fSaveEmerging = false;


  //PrimaryGenerator Neutron energy output
  Double_t fbranchEInitN {-1};
  //Event Number
  Int_t fEventNo {-1};
  //EventAction scored
  Bool_t fbranchValid = false;
  //secondaryParticles
  std::vector<Double_t> fbranchSecondaryNeutron;
  std::vector<Double_t> fbranchSecondaryNeutronEmerging;
  std::vector<Double_t> fbranchSecondaryElectron;
  std::vector<Double_t> fbranchSecondaryElectronEmerging;
  std::vector<Double_t> fbranchSecondaryGamma;
  std::vector<Double_t> fbranchSecondaryGammaEmerging;
  //FissionProducts
  std::vector<std::string> fbranchFissionProduct;
  std::vector<Double_t> fbranchFissionProductMass;
  std::vector<Int_t> fbranchFissIon;
  std::vector<Double_t> fbranchFissNeut;
  std::vector<Double_t> fbranchFissNeutEmerging;


  // Emerging particles
    // Emerging particle record
  std::vector<Int_t> fbranchEmergingRank;
  std::vector<Int_t> fbranchEmergingThreadId;
  std::vector<Int_t> fbranchEmergingEventId;
  std::vector<Int_t> fbranchEmergingId;
  std::vector<Int_t> fbranchEmergingParentId;
  std::vector<Int_t> fbranchEmergingPDG;
  std::vector<ROOT::Math::XYZTVector> fbranchEmergingPos;
  std::vector<ROOT::Math::XYZTVector> fbranchEmergingP; 
  std::vector<std::string> fbranchEmergingProcess; 


  //Neutron Energy spectrum
  std::vector<Double_t> fEmergingNeutrons;

  std::vector<Double_t> fbranchNEmissionSpec;
  std::vector<ROOT::Math::XYZVectorD> fbranchNEmissionSpecVer;

  //Trees parameters
  G4String FileOutName = "-Spectrum";
  G4String FileOutNameRank;
  G4int fThreadid {-10};
  std::vector <TFile*> fout {0};
  std::vector <TTree*> tree {0};

  G4String EmergingFileOutName = "-EmergingParticles";
  G4String EmergingFileOutNameRank;
  std::vector <TFile*> fEmergingOut {0};
  std::vector <TTree*> fEmergingTree {0};

  //AmBe Cross section data
  TGraph *gXS_0;
  TGraph *gXS_1;
  TGraph *gXS_2;
  TGraph *gXS_3;
  TGraph *gstopping;
  TGraph *gXS_t;
  //AmBe Cross section loader
  TGraph2DErrors* loadAng(TString filename);
  void loadFile(TString filename, std::vector<double>& Z_array, std::vector<double>& X_array, std::vector<double>& Y_array, std::vector<double>& ex_array, std::vector<double>& ey_array);
  std::vector<TGraph*> loadLegendre(TString filename);
  //TGraph2DErrors *g_gs;
  //TGraph2DErrors *gp;
  //TGraph2DErrors *gpp;
  std::vector<TGraph*> g_gs;
  std::vector<TGraph*> gp;
  std::vector<TGraph*> gpp;
  std::vector<TGraph*> gppp;
  Int_t legendrePolLimit = 8;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
