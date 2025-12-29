#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "NNDCLoader.hh"

#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <G4Geantino.hh>
#include <G4ParticleGun.hh>
#include <Randomize.hh>
#include <G4GeneralParticleSource.hh>
#include <G4Threading.hh>
#include <TROOT.h>
#include "G4AutoLock.hh"

#include <TMath.h>
#include <TF1.h>

namespace
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction, DetectorConstruction*  detectorConstruction, G4String IsotopeString, G4int NumberOfThreads)
: G4VUserPrimaryGeneratorAction(),
  fGPS(0),
  fEnvelopeBox(0),
  fRunAction(runAction),
  fDetectorConstruction(detectorConstruction),
  fRadioIsotope(IsotopeString),
  fNumberOfThreads(NumberOfThreads)
{
  ROOT::EnableThreadSafety();//see https://root-forum.cern.ch/t/mutexes-when-running-inside-geant4-threads/42060
  G4AutoLock l(&aMutex);

  //messenger
  //f241AmGammaEnabled = true;
  f241AmGammaEnabled = false;
  PrimaryGeneratorMessenger();

  //G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  //G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("electron");
  G4ParticleDefinition* geantino = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  G4ParticleDefinition* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");

  //GPS and ThreadSafe Init
  fGPS.resize(fNumberOfThreads, nullptr);
  fNeutronGPS.resize(fNumberOfThreads, nullptr);
  fCarbonGPS.resize(fNumberOfThreads, nullptr);
  fGamma_13_9.resize(fNumberOfThreads, nullptr);
  fGamma_26_34.resize(fNumberOfThreads, nullptr);
  fGamma_59_54.resize(fNumberOfThreads, nullptr);
  
  fThreadid = G4Threading::G4GetThreadId();//omp_get_thread_num();
  G4cout << "PrimaryGenerator:: ThreadID "<< fThreadid << G4endl;
  if(fThreadid<0) fThreadid=0;
  if (fThreadid<0 or fThreadid>fNumberOfThreads-1) G4cout << "ERROR AT THREAD NUMBER: " << fThreadid << G4endl;
  
  if (true)
  {
    fGPS[fThreadid] = new G4GeneralParticleSource();
     // Defensive: keep only source 0
    {
      int nsrc = fGPS[fThreadid]->GetNumberofSource();
      for (int k = nsrc - 1; k >= 1; --k) fGPS[fThreadid]->DeleteaSource(k);
      fGPS[fThreadid]->SetCurrentSourceto(0);
    }

    fGPS[fThreadid]->GetCurrentSource()->SetNumberOfParticles(1);
    fGPS[fThreadid]->GetCurrentSource()->SetParticleDefinition(neutron);
    fNeutronGPS[fThreadid] = fGPS[fThreadid]->GetCurrentSource();
    fNeutronGPS[fThreadid]->GetPosDist()->SetPosDisType("Point");
    fNeutronGPS[fThreadid]->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
    //fNeutronGPS[fThreadid]->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));//DEBUG
    fNeutronGPS[fThreadid]->GetAngDist()->SetAngDistType("iso");
    fNeutronGPS[fThreadid]->GetEneDist()->SetEnergyDisType("Mono");
    fGPS[fThreadid]->SetMultipleVertex(false);

    
    if (fRadioIsotope!="USF" and (fRadioIsotope.length() < 6 or fRadioIsotope.substr(0,6)!="Single"))
    {
      fGPS[fThreadid]->AddaSource(1);
      fGPS[fThreadid]->GetCurrentSource()->SetNumberOfParticles(1);
      fCarbonGPS[fThreadid] = fGPS[fThreadid]->GetCurrentSource();
      fCarbonGPS[fThreadid]->SetParticleDefinition(geantino);
      fCarbonGPS[fThreadid]->GetPosDist()->SetPosDisType("Point");
      fCarbonGPS[fThreadid]->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
      //fCarbonGPS[fThreadid]->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));//DEBUG
      fCarbonGPS[fThreadid]->GetAngDist()->SetAngDistType("iso");
      fCarbonGPS[fThreadid]->GetEneDist()->SetEnergyDisType("Mono");

      if (f241AmGammaEnabled)
      {
        if (alphaPerNeutron<=0)
        {
          G4cerr << "PrimaryGenerator::Init alphaPerNeutron not defined"<<G4endl;
          return;
        }
        fGPS[fThreadid]->AddaSource(1);
        fGamma_59_54[fThreadid] = fGPS[fThreadid]->GetCurrentSource();
        fGPS[fThreadid]->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.37);
        fGamma_59_54[fThreadid]->SetParticleDefinition(gamma);
        fGamma_59_54[fThreadid]->GetPosDist()->SetPosDisType("Point");
        fGamma_59_54[fThreadid]->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
        fGamma_59_54[fThreadid]->GetAngDist()->SetAngDistType("iso");
        fGamma_59_54[fThreadid]->GetEneDist()->SetEnergyDisType("Mono");
        fGamma_59_54[fThreadid]->GetEneDist()->SetMonoEnergy(59.5409*keV);

        fGPS[fThreadid]->AddaSource(1);
        fGamma_26_34[fThreadid] = fGPS[fThreadid]->GetCurrentSource();
        fGPS[fThreadid]->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.0227);
        fGamma_26_34[fThreadid]->SetParticleDefinition(gamma);
        fGamma_26_34[fThreadid]->GetPosDist()->SetPosDisType("Point");
        fGamma_26_34[fThreadid]->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
        fGamma_26_34[fThreadid]->GetAngDist()->SetAngDistType("iso");
        fGamma_26_34[fThreadid]->GetEneDist()->SetEnergyDisType("Mono");
        fGamma_26_34[fThreadid]->GetEneDist()->SetMonoEnergy(26.3446*keV);

        fGPS[fThreadid]->AddaSource(1);
        fGamma_13_9[fThreadid] = fGPS[fThreadid]->GetCurrentSource();
        fGPS[fThreadid]->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.359);
        fGamma_13_9[fThreadid]->SetParticleDefinition(gamma);
        fGamma_13_9[fThreadid]->GetPosDist()->SetPosDisType("Point");
        fGamma_13_9[fThreadid]->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
        fGamma_13_9[fThreadid]->GetAngDist()->SetAngDistType("iso");
        fGamma_13_9[fThreadid]->GetEneDist()->SetEnergyDisType("Mono");
        fGamma_13_9[fThreadid]->GetEneDist()->SetMonoEnergy(13.9*keV);
      }
    }
    fGPS[fThreadid]->SetMultipleVertex(true);
  }
  else
  {
    // No carbon/gammas in USF/Single
    fCarbonGPS[fThreadid] = nullptr;
    fGamma_59_54[fThreadid] = nullptr;
    fGamma_26_34[fThreadid] = nullptr;
    fGamma_13_9[fThreadid]  = nullptr;
  }

  if (fRadioIsotope!="USF" and (fRadioIsotope.length() < 6 or fRadioIsotope.substr(0,6)!="Single"))
  {
    //Hard-code diff Xsection for 9Be(a,n)12C - Part of PrimaryGenerator
    //renaming these is useless as this simply gets the pointer to something else
    if (true)
    {
      gXS_0 = fRunAction->GetgXS_0();
      gXS_1 = fRunAction->GetgXS_1();
      gXS_2 = fRunAction->GetgXS_2();
      gXS_3 = fRunAction->GetgXS_3();
      gstopping = fRunAction->Getgstopping();
      gXS_t = fRunAction->GetgXS_t();
      g_gs = fRunAction->Getg_gs();
      gp = fRunAction->Getgp();
      gpp = fRunAction->Getgpp();
      gppp = fRunAction->Getgppp();

      fIntXSMax = 1.2*(TMath::MaxElement(gXS_0->GetN(),gXS_0->GetY())+TMath::MaxElement(gXS_1->GetN(),gXS_1->GetY())+TMath::MaxElement(gXS_2->GetN(),gXS_2->GetY()));
      G4cout << "PrimaryGeneratorAction::Maximum integrated XS*1.2 = " <<fIntXSMax<<G4endl;
      fLeg8 = new TF1("legendre8","sin(x)*([0]+[1]*cos(x)+[2]*(3*pow(cos(x),2)-1)/2+[3]*(5*pow(cos(x),3)-3*(cos(x)))/2+[4]*(35*pow(cos(x),4)-30*pow(cos(x),2)+3)/8+[5]*(63*pow(cos(x),5)-70*pow(cos(x),3)+15*(cos(x)))/8+[6]*(231*pow(cos(x),6)-315*pow(cos(x),4)+105*pow(cos(x),2)-5)/16+[7]*(429*pow(cos(x),7)-693*pow(cos(x),5)+315*pow(cos(x),3)-35*(cos(x)))/16+[8]*(6435*pow(cos(x),8)-12012*pow(cos(x),6)+6930*pow(cos(x),4)-1260*pow(cos(x),2)+35)/128)",0,TMath::Pi());
    }

    if (fRadioIsotope=="241Am")
      fNNDCAlpha = new NNDCLoader("AmBeData/241AmAlpha.dat");
    else if (fRadioIsotope=="239Pu")
      fNNDCAlpha = new NNDCLoader("AmBeData/239PuAlpha.dat");

    r_max = fDetectorConstruction->GetAmBeSolid()->GetOuterRadius() - max_depth*mm;
    z_max = fDetectorConstruction->GetAmBeSolid()->GetZHalfLength() - max_depth*mm;
  }
  else if (fRadioIsotope.length()>=6 and fRadioIsotope.substr(0,6)=="Single")
  {
    //singleNeutronEnergy = std::stod(fRadioIsotope.substr(6,fRadioIsotope.length()))*MeV;
    try
    {
      singleNeutronEnergy = std::stod(fRadioIsotope.substr(6)) * MeV;
    }
    catch(const std::exception& e)
    {
      G4cerr << "PrimaryGeneratorAction::Init Bad Single energy string: " << fRadioIsotope << " : " << e.what() << G4endl;
      singleNeutronEnergy = -1;
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if (fGPS[fThreadid]) delete fGPS[fThreadid];
  if (fNNDCAlpha) delete fNNDCAlpha;
  if (fMessenger) delete fMessenger;
  if (fLeg8) delete fLeg8;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fThreadid = G4Threading::G4GetThreadId();//omp_get_thread_num();
  if(fThreadid<0) fThreadid=0;
  // if (fThreadid<0 || fThreadid>fNumberOfThreads-1) G4cout << "ERROR AT THREAD NUMBER: " << fThreadid << G4endl;

  if (fRadioIsotope=="USF")
  {
    G4double EnerU = USF();
    fNeutronGPS[fThreadid]->GetEneDist()->SetMonoEnergy(EnerU);
    //Save initital Neutron energy in the Ntuple
    fRunAction->FillInitialNeutron(EnerU);
    fGPS[fThreadid]->GeneratePrimaryVertex(anEvent);
    return;
  }
  if (fRadioIsotope.length()>=6 and fRadioIsotope.substr(0,6)=="Single")
  //else if (singleNeutronEnergy!=-1)
  {
    fNeutronGPS[fThreadid]->GetEneDist()->SetMonoEnergy(singleNeutronEnergy);
    //Save initital Neutron energy in the Ntuple
    fRunAction->FillInitialNeutron(singleNeutronEnergy);
    fGPS[fThreadid]->GeneratePrimaryVertex(anEvent);
    return;
  }
  else
  {
  
    //Randomise Centre
    G4double z = (2*G4UniformRand()-1) * z_max;
    G4double r = std::sqrt(G4UniformRand()) * r_max;
    G4double theta = G4UniformRand()*2*pi;

    G4double x, y;
    x = cos(theta)*r;
    y = sin(theta)*r;
    G4ThreeVector PositionVector = G4ThreeVector(x,y,z);

    fNeutronGPS[fThreadid]->GetPosDist()->SetCentreCoords(PositionVector);
    fCarbonGPS[fThreadid]->GetPosDist()->SetCentreCoords(PositionVector);

    if (f241AmGammaEnabled)
    {
      G4ThreeVector GammaVector = PositionVector;//G4ThreeVector(cos(theta)*r_max,sin(theta)*r_max,z);
      fGamma_59_54[fThreadid]->GetPosDist()->SetCentreCoords(GammaVector);
      fGamma_26_34[fThreadid]->GetPosDist()->SetCentreCoords(GammaVector);
      fGamma_13_9[fThreadid]->GetPosDist()->SetCentreCoords(GammaVector);
    }

    /* 241Am
    * 84.8% 5.48556 MeV
    * 13.1% 5.44280 MeV
    * Tot: 97.9%. Re-averaged: 87.62% and 13.38%
    */
    Double_t Einitial = 5.48556;//Main peak
    //if(G4UniformRand()<0.1338) Einitial = 5.44280;//Satellite peak
    Einitial = fNNDCAlpha->GetEnergy(G4UniformRand());
    if (fDBG) G4cout << "PrimaryGenerator: Einitial="<<Einitial<<G4endl;//Debug
    /* //DBG - Check ratios
    if (true)
    {
      std::ofstream NNDCFile("NNDCOutput.dat", std::ios::out | std::ios::app);
      if (NNDCFile)
      {
        NNDCFile << Einitial<<"\n";
      }
      NNDCFile.close();
    }
    //*/
    if (Einitial==-1)
    {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "NNDC001", EventMustBeAborted, "Invalid energy provided from NNDCLoader");
      return;
    }
    Double_t alphaEnergy = InteractionE(Einitial);//alpha of energy when it interacts - after Eloss

    G4double EnerNeutron = (GetNeutron(alphaEnergy))*MeV;
    if (EnerNeutron<0.)
    {
      G4cerr << "PrimaryGenerator: Negative neutron energy: "<<EnerNeutron<<G4endl;
      return;
    }
    fNeutronGPS[fThreadid]->GetEneDist()->SetMonoEnergy(EnerNeutron);
    
//    if (true)//Test single states only
    if (false)
    {
      //if (excited or excited2 or excited3)
      //if (!excited or excited2 or excited3)
      //if (excited or !excited2 or excited3)
      //if (excited or excited2 or !excited3)
      //return;
    }

    //Save initital Neutron energy in the Ntuple
    fRunAction->FillInitialNeutron(EnerNeutron);
    //do not close row toa void data loss (otherwise try with closing every row, i.e. not conditioned on activation of detectors in EventAction)

    Double_t _CarbonEnergy = alphaEnergy+Q_9BeAN-EnerNeutron/MeV;//-ExcitationEnergy
    Double_t CarbonEnergy = _CarbonEnergy*MeV;
    if (CarbonEnergy<0.*keV)
    {
      if (CarbonEnergy>-1.0*keV)//0>CarbonEnergy>-1.0keV - acceptable range, set to 0
      {CarbonEnergy=0.0*keV;}
      else
      {
        G4cout << "PrimaryGenerator: ";
        if (excited) G4cout << "Carbon in 1st excited state; ";
        if (excited2) G4cout << "Carbon in 2nd excited state; ";
        else G4cout << "Carbon in ground state; ";
        G4cerr << "Carbon Energy negative: "<< CarbonEnergy << "; neutron energy:"<< EnerNeutron <<G4endl;
      }
    }
    //Gamma from Carbon:
    //First test with just the gamma, then create a proper carbon (because of doppler broadening)
    //Also, for 2nd excited, check if hoyle state decays by alpha or just gamma. If just gamma, apply the branching ratio to gamma by creating this state using G4UniformRand
    if (excited or excited2 or excited3)
    {
      //fGPS[fThreadid]->SetMultipleVertex(false);
      if (excited)
      {
        if (CarbonEnergy-C_Excit_1ST<0) G4cerr << "PrimaryGenerator: Carbon energy below 1st excited"<<G4endl;
        //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(4.43982*MeV);
        G4ParticleDefinition* C_Ion_1ST = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_1ST);
        fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_1ST);
        fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
        fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_1ST);
        if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy-C_Excit_1ST)/MeV<<" MeV, 1st excited"<< G4endl;//Debug
      }
      else if (excited2)
      {
        if (CarbonEnergy-C_Excit_2ND<0) G4cerr << "PrimaryGenerator: Carbon energy below 2nd excited"<<G4endl;
        if(G4UniformRand()<4.16E-4)//gamma from Hoyle
        //if (true)
        {
          //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(7.65407*MeV);
          G4ParticleDefinition* C_Ion_2ND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_2ND);
          fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_2ND);
          fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
          fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_2ND);
          if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy-C_Excit_2ND)/MeV<<" MeV, 2nd excited"<< G4endl;//Debug
        }
        ///*
        else//alpha decay
        {
          //set as null carbon as it should alpha decay
          G4ParticleDefinition* C_Ion_GND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_GND);
          fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_GND);
          fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
          fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_2ND);
          //G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy)/MeV<<" MeV, GND state"<< G4endl;//Debug
          //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(0);
          if (fDBG) G4cout << "PrimaryGenerator: Null Carbon 3RD"<< G4endl;//Debug
        }
        //*/
      }
      else if (excited3)
      {
        if (CarbonEnergy-C_Excit_3RD<0) G4cerr << "PrimaryGenerator: Carbon energy below 3rd excited"<<G4endl;
        if(G4UniformRand()<4.1E-7)//gamma
        //if (true)
        {
          //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(7.65407*MeV);
          G4ParticleDefinition* C_Ion_3RD = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_3RD);
          fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_3RD);
          fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
          fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_3RD);
          if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy-C_Excit_3RD)/MeV<<" MeV, 3rd excited"<< G4endl;//Debug
        }
        ///*
        else//alpha decay
        {
          //set as null carbon as it should alpha decay
          G4ParticleDefinition* C_Ion_GND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_GND);
          fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_GND);
          fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
          fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_3RD);
          //G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy)/MeV<<" MeV, GND state"<< G4endl;//Debug
          //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(0);
          if (fDBG) G4cout << "PrimaryGenerator: Null Carbon 3RD"<< G4endl;//Debug
        }
        //*/
      }
      else
      {
        G4cerr << "PrimaryGenerator: ERROR! - flags have an issue in PrimaryGeneratorAction::GeneratePrimaries()"<<G4endl;
        return;
      }
      //if (excited) G4cout << "Excited 1st , Primaries : " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
      //else if (excited2) G4cout << "Excited 2nd , Primaries : " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
    }
    else //GND state
    {
      //fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(0.0*MeV);
      G4ParticleDefinition* C_Ion_GND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_GND);
      fCarbonGPS[fThreadid]->SetParticleDefinition(C_Ion_GND);
      fCarbonGPS[fThreadid]->SetParticleCharge(ionCharge);
      fCarbonGPS[fThreadid]->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_GND);
      if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy)/MeV<<" MeV, GND state"<< G4endl;//Debug
    }
    fGPS[fThreadid]->GeneratePrimaryVertex(anEvent);
    return;
    //G4cout << "Primaries: " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
    //if (anEvent->GetNumberOfPrimaryVertex()!=2) G4cout << "ERROR!" << G4endl; //-- All good //DBG
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::GetNeutron(Double_t alphaEnergy)
{
  Double_t Eint = alphaEnergy;//energy of interacting alpha
  //Sample theta CM
  Double_t thetaCM = InteractionTheta(Eint);
  if(thetaCM==0) return 0;//Skip this guy?
  else if (thetaCM == -1)
  {
    G4cout << "PrimaryGenerator::GetNeutron\tError when determining the reaction angle in the Centre of Mass frame. Sampling failed in InteractionTheta"<<G4endl;
    return -1;
  }
  //thetalab = ConvertCMtoLab(Eint, thetaCM);
  Double_t En = NeutronEnergy(Eint, thetaCM);
  if(En==0) return 0;//Skip this guy
  return En;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PrimaryGeneratorAction::InteractionE(Double_t Eint)
{
  //Rejection sampling for z_interaction based on XS(Eint) and stopping power
  //Then return Eint for z
  Bool_t trials=true;
  //Double_t XSmax = 720;
  Double_t XSmax = fIntXSMax;
  //G4cout << "XSMax " <<XSmax<<G4endl;
  /*
   * XSmax  x1    x1.2
   * 0+1+2  573   687
   * ...+3  638   766
   * ...+b  1072  1286
   */
  Int_t loopcounter=0;
  while(trials) //Randomly choose a depth into the target (rather than energy to keep same target thickness)
  {
    Double_t trial_z = max_depth*G4UniformRand();
    //Calc XS at Z
    Double_t energy = Eloss(Eint,trial_z);
    if(energy<=Eth) continue;
    //Double_t XS = (gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
    //Double_t XS = (gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy)+gXS_3->Eval(energy));
    Double_t XS = gXS_t->Eval(energy);
    Double_t choose_chan = G4UniformRand();
    Double_t BRn0=gXS_0->Eval(energy)/XS;
    Double_t BRn1=gXS_1->Eval(energy)/XS;
    Double_t BRn2=gXS_2->Eval(energy)/XS;
    Double_t BRn3=gXS_3->Eval(energy)/XS;
    excited=false;
    excited2=false;
    excited3=false;

    //G4cout << "PrimaryGeneratorAction::InteractionE - select state: "<<Eint<<" "<<choose_chan<<" "<<BRn0<<" "<<BRn1<<" "<<BRn2<<" "<<BRn3<<" "<<XS<<G4endl;

    if(choose_chan<BRn0)
    {
      excited=false;
      excited2=false;
      excited3=false;
    }
    else if(choose_chan<BRn0+BRn1)
    {
      excited=true;
      excited2=false;
      excited3=false;
    }
    else if(choose_chan<BRn0+BRn1+BRn2)
    {
      excited=false;
      excited2=true;
      excited3=false;
    }
    else
    {
      excited=false;
      excited2=false;
      excited3=true;
    }

    if(XS>XSmax) {G4cout<<"PrimaryGenerator::IneractionE\tMax XS is exceeded!!"<<G4endl;}
    if(XS>XSmax*G4UniformRand())
    {
      return energy;
    }
    loopcounter++;
    if(loopcounter>1000)
    {
      G4cerr<<"PrimaryGenerator::InteractionE\tRejection sampling not working properly - 1000 samples taken"<<G4endl;
      return -1;
    }
  }
  G4cerr<<"PrimaryGenerator::InteractionE\tRejection sampling exited - this should not occour"<<G4endl;
  return -2;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PrimaryGeneratorAction::Eloss(Double_t Ea, Double_t z) //Energy of alpha of initial energy Ea after z mm in Be
{
  /*
   * Manually extracted from Lise++
   * Using 241Am at stechio 110
   * 16O at stechio 220
   * 9Be at stechio 890
   *
   * and alpha energy at max-> 5.48556
   *
   * -> max range 36.82486*um;
  */
  Double_t Eout, range;
  if (fRadioIsotope=="241Am")
  {
    //Double_t range=0.000528638*Ea*Ea+0.00355149*Ea+0.000974915;//mm
    range=0.000530302*Ea*Ea+0.00354792*Ea+0.00105601;//mm
    if(z>range) return 0;//Fully stopped before Z
    range-=z;
    //Double_t Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
    Eout=-1591.4*range*range+211.143*range-0.0804195;//MeV
  }
  else if (fRadioIsotope=="239Pu")
  {
    //Double_t range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
    range=0.000530474*Ea*Ea+0.00355288*Ea+0.00105571;//mm
    if(z>range) return 0;//Fully stopped before Z
    range-=z;
    //Double_t Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
    Eout=-1609.33*range*range+211.354*range-0.0815178;//MeV
  }
  if (Eout<0.0)
  {
    //G4cerr << "PrimaryGenerator::Eloss yields negative energy: "<<Eout<<G4endl;
    Eout = 0.;
  }
  return Eout;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PrimaryGeneratorAction::InteractionTheta(Double_t E_inc) //COM theta
{
  Double_t random_theta = acos(-1+2*G4UniformRand());
  //random_theta = 2.*pi*G4UniformRand();//NOTE:check this
  //Double_t random_theta = pi*G4UniformRand();//TEST -> It is effetively as in two lines before acos(-1+2*G4UniformRand());

  Bool_t sample=true;
  Double_t XSmax=1.2;//Max is 1, XS data are in realtive yield
  Int_t counter=0;
  Double_t A[9] = {0.};//initialise to 0
  Double_t diffXSNormalise = 1;
  std::vector<TGraph*> select;

  if (!excited and !excited2 and !excited3) select=g_gs;
  else if (excited and !excited2 and !excited3) select=gp;
  else if (!excited and excited2 and !excited3) select=gpp;
  else if (!excited and !excited2 and excited3) select=gppp;
  else G4cerr << "PrimaryGenerator: ERROR! - flags have an issue in PrimaryGeneratorAction::InteractionTheta()"<<G4endl;

  //G4cout << "PrimaryGeneratorAction::InteractionTheta - loading Legendre pol for "<<excited<<excited2<<excited3<<" with Ea="<<E_inc<<G4endl;

  for (int i=0; i<=8; i++)
  {
    if (select[i])
    {
      //G4cout <<"Exists "<<i<<G4endl;
      Double_t val = select[i]->Eval(E_inc);
      A[i] = val;
      fLeg8->SetParameter(i,val);
    }
    else
    {
      G4cerr << "LegendrePol is nullptr"<<G4endl;
      A[i] = 0;
    }
  }

  diffXSNormalise = fLeg8->GetMaximum();
  if (diffXSNormalise <= 0.0)
  {
    G4cerr << "PrimaryGenerator: diffXSNormaliseN0 invalid, defaulting to 1" << G4endl;
    diffXSNormalise = 1.0;
  }
  if (fDBG>3) G4cout << "LegPols\n"<<A[0]<<" "<<A[1]<<" "<<A[2]<<" "<<A[3]<<" "<<A[4]<<" "<<A[5]<<" "<<A[6]<<" "<<A[7]<<" "<<A[8]<<" "<<G4endl;

  while(sample)
  {
    Double_t XS_sample = XSmax*G4UniformRand();
    Double_t XS = -1;

    //Use this to use Legendre Polynomials as in GVDZ 1976(?) as reported on exfor https://www-nds.iaea.org/exfor//servlet/X4sGetSubent?subID=130001004
    random_theta = 2.*pi*G4UniformRand();
    Double_t cos_RNDtheta = cos(random_theta);
    XS = (sin(random_theta)*(A[0]+A[1]*cos_RNDtheta+A[2]*(3*pow(cos_RNDtheta,2)-1)/2+A[3]*(5*pow(cos_RNDtheta,3)-3*(cos_RNDtheta))/2+A[4]*(35*pow(cos_RNDtheta,4)-30*pow(cos_RNDtheta,2)+3)/8+A[5]*(63*pow(cos_RNDtheta,5)-70*pow(cos_RNDtheta,3)+15*(cos_RNDtheta))/8+A[6]*(231*pow(cos_RNDtheta,6)-315*pow(cos_RNDtheta,4)+105*pow(cos_RNDtheta,2)-5)/16+A[7]*(429*pow(cos_RNDtheta,7)-693*pow(cos_RNDtheta,5)+315*pow(cos_RNDtheta,3)-35*(cos_RNDtheta))/16+A[8]*(6435*pow(cos_RNDtheta,8)-12012*pow(cos_RNDtheta,6)+6930*pow(cos_RNDtheta,4)-1260*pow(cos_RNDtheta,2)+35)/128))/diffXSNormalise;//normalise by maximum of diffXS for the specified energy between 0 and pi

    if (XS==-1) G4cerr << "ERROR XS not set: -1"<<G4endl;
    //*/

    if (XS>XSmax) G4cerr<<"PrimaryGenerator::InteractionTheta\tSampling error for angle: XS "<<XS<<" bigger than XSmax "<<XSmax<<". Ealpha = "<<E_inc<<"\t"<<excited<<"\t"<<excited2<<"\t"<<excited3<<G4endl;
    //if (XS==0 or excited2) XS=1;//Outside of the energy range for this data set or Hoyle (poor data)
    if (XS==0)
    {
      G4cout << "PrimaryGenerator: XS is zero. Invalid, defaulting to 1" << G4endl;
      XS=1;//Outside of the energy range for this data set (poor data)
    }
    if (XS_sample<XS)
    {
      return random_theta;//Isotropic for now
    }
    counter++;
    if (counter>100)
    {
      G4cerr<<"PrimaryGenerator::InteractionTheta\tOver 100 samples taken - no solution for XS "<<XS<<" and XS_Sample "<<XS_sample<<", Ealpha "<<E_inc<<"\t"<<excited<<"\t"<<excited2<<"\t"<<excited3<<G4endl;
      return 0;
    }
  }
  //Should not reach here
  G4cerr << "PrimaryGenerator::InteractionTheta\tError in determination - generic error"<<G4endl;
  return -1;//should not reach //TODO:check if correct
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PrimaryGeneratorAction::ConvertCMtoLab(Double_t Eint, Double_t thetaCM) //Convert thetaCM to lab
{
  Double_t Ex=C_Excit_GND_MeV;
  if(excited) Ex=C_Excit_1ST_MeV;
  if(excited2) Ex=C_Excit_2ND_MeV;
  if(excited3) Ex=C_Excit_3RD_MeV;

  Double_t gamma = sqrt(((m_4He*m_n)/(m_9Be*m_12C))*(Eint/(Eint+(Q_9BeAN-Ex)*(1.+m_4He/m_9Be))));//Calculate gamma for the conversion to lab COM
  gamma = gamma*((m_4He+m_9Be)/(m_n+m_12C));//Account for change in the CM velocity
  Double_t thetalab = atan2(sin(thetaCM),(cos(thetaCM)+gamma));//
  secondks = false;
  if(gamma*cos(thetaCM)<=-1.) secondks=true;//Neutron is going bac
  //thetalab=(12*G4UniformRand())*d2r;//TEST
  //if(G4UniformRand()>0.5) secondks=true;//TEST
  return thetalab;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PrimaryGeneratorAction::NeutronEnergy(Double_t Eb, Double_t thetaCM) //Get neutron energy for a given Eint and theta
{
  Double_t En;
  Double_t Ex=C_Excit_GND_MeV;
  if(excited) Ex=C_Excit_1ST_MeV;//MeV
  if(excited2) Ex=C_Excit_2ND_MeV;
  if(excited3) Ex=C_Excit_3RD_MeV;

  Double_t ECM = m_9Be*Eb/(m_9Be+m_4He)+Q_9BeAN-Ex;
  if(ECM<0)
  {
    G4cerr<<"PrimaryGenerator::NeutronEnergy\tBelow threshold! Something is wrong!"<<G4endl;
    return 0;
  }
  Double_t p_nT = sqrt(2.*ECM*m_n*(m_12C/(m_n+m_12C)));
  Double_t phi=2*pi*G4UniformRand();
  Double_t pn[3]={p_nT*sin(thetaCM)*cos(phi),p_nT*sin(thetaCM)*sin(phi),p_nT*cos(thetaCM)};//COM neutron momentum
  pn[2]+=(m_n/(m_12C+m_n))*sqrt(2.*m_4He*Eb);//Add the momentum boost for the lab frame
  Double_t En_new = 0;
  for(int i=0;i<3;i++) En_new+=pn[i]*pn[i];
  En_new /= (2.*m_n);
  En = En_new;
  return En;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::PrimaryGeneratorMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/primaryGenerator/", "primaryGenerator control");
  fMessenger->DeclareProperty("241AmGamma",f241AmGammaEnabled);
  fMessenger->DeclareProperty("radioIsotope",fRadioIsotope);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
