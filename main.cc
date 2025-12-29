#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "PhysicsList.hh"

#include <vector>
#include <mpi.h>
#include <unistd.h>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include <G4RunManagerFactory.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>
#include <G4String.hh>
#include <G4UImanager.hh>

#include <QGSP_BERT_HP.hh>
#include <QGSP_BIC_HP.hh>
#include <QGSP_BIC_HPT.hh>
#include <G4ParticleHPManager.hh>
#include <G4PhysListFactory.hh>
#include <G4ParallelWorldPhysics.hh>

#include <G4MPImanager.hh>
#include <G4MPIsession.hh>


#include <G4HadronPhysicsShielding.hh>
#include <Randomize.hh>

G4long setUniqueSeed(int rank);

int main(int argc, char** argv)
{
  setenv("QT_STYLE_OVERRIDE", "Fusion", 1);//Avoid Qt related segfault when using Breeze theme

  G4cout << "Application starting..." << G4endl;

  G4bool MPI_ = false;
  if (const char* env_mpi = std::getenv("AMBE_MPI"))
  {
    G4String mpi_str(env_mpi);
    if (mpi_str=="1" or mpi_str=="true" or mpi_str=="TRUE" or mpi_str=="True")
    {
      MPI_=true;
    }
    else if (mpi_str=="0" or mpi_str=="false" or mpi_str=="FALSE" or mpi_str=="False")
    {
      MPI_=false;
    }
    else
    {
      MPI_=false;
    }
  }
  G4cout << "Use MPI: "<<MPI_<<G4endl;

  std::vector<G4String> macros;
  G4bool interactive = false;

  G4int NumberOfThreads = 1;
  G4int noWaterBath = 0;
  /*
   * 0 is water bath
   * 1 is no water
   * 2 is world is of Water
   */
  G4bool FissFragments=false;
  G4bool ScoreGamma=false;
  G4int Isotope = -1;
  G4String IsotopeString;
  /*
   * Isotope selector:
   * -1 Default
   * 0 AmBe
   * 1 ^{239}PuBe
   * 2 SF for U
   * 3 MonoEnergetic Neutrons
   */
  G4bool InitialNeutrons = false;
  G4int CasingSelection = 0;
  /*
   * Casing selectot:
   * 0 cylindrical approximation
   * 1 X3 casing
   */
  G4bool AzimuthalScoring = false;


  // random engine
  CLHEP::MTwistEngine randomEngine;
  G4Random::setTheEngine(&randomEngine);

  std::vector<G4String> arguments;

  G4MPImanager* g4MPI = nullptr;
  if (MPI_)
  {
    g4MPI = new G4MPImanager(argc,argv);
    arguments = g4MPI->ReturnArguments();
  }
  else
  {
    for (int i = 1; i < argc; i++)
    {
      arguments.push_back(argv[i]);
      //G4cout << "Argument " << i << ": " << arguments[i-1] << G4endl;

    }
  }

  if (!MPI_ and argc==1)
  {
    interactive = true;
  }
  else
  {
    for (unsigned int i = 0; i < arguments.size(); i++)
    {
      G4String arg = arguments[i];
      if (arg == "-i" || arg == "--interactive")
      {
        interactive = true;
      }
      else if (arg == "-t" || arg == "--threads")
      {
        NumberOfThreads = std::stoi(arguments[i+1]);
        i++;
      }
      else if (arg == "-nw" || arg == "--nowater")
      {
        noWaterBath = std::stoi(arguments[i+1]);
        i++;
      }
      else if (arg == "-ffs" || arg == "--fissionFragmentsScore")
      {
        FissFragments=true;
      }
      else if (arg == "-r" || arg == "--isotope")
      {
        if (arguments[i+1].length()>=6 and arguments[i+1].substr(0,6)=="Single")
        {
          Isotope = 3;
          IsotopeString = arguments[i+1];
        }
        else
        {
          Isotope = std::stoi(arguments[i+1]);
        }
        i++;
      }
      else if (arg == "-in" || arg == "--initialNeutrons")
      {
        InitialNeutrons = true;
      }
      else if (arg == "-cs" || arg == "--casing")
      {
        CasingSelection = std::stoi(arguments[i+1]);
        i++;
      }
      else if (arg == "-ass" || arg == "--azimuthalSurfaceScoring")
      {
        AzimuthalScoring = true;
      }
      else
      {
        if (!MPI_)
        {
          macros.push_back(arg);
        }
      }
    }
  }
  if (Isotope<=0) IsotopeString = "241Am";
  else if (Isotope==1) IsotopeString = "239Pu";
  else if (Isotope==2) IsotopeString = "USF";


  // MPI session (G4MPIsession) instead of G4UIterminal
  G4MPIsession* g4MPISession  = nullptr;
  if (MPI_)
  {
    g4MPISession = g4MPI->GetMPIsession();
    // LAM/MPI users can use G4tcsh.
    G4String prompt = "[40;01;33m";
    prompt += "G4MPI";
    prompt += "[40;31m(%s)[40;36m[%/][00;30m:";
    g4MPISession->SetPrompt(prompt);
  }

  // Create the run manager
  //G4cout <<"Initialising threading" <<G4endl;
  auto* runManager  =
    //G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    //G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
    //G4RunManagerFactory::CreateRunManager(G4RunManagerType::TBB);//BEAR
  //G4MTRunManager* runManager = new G4MTRunManager();
  runManager->SetVerboseLevel(1);
  runManager->SetNumberOfThreads(NumberOfThreads);
  G4cout << "THREADS COUNT " << runManager->GetNumberOfThreads() << "\n";

  G4int actualNumberOfThreads = runManager->GetNumberOfThreads();
  /*
  G4int actualNumberOfThreads = -1;
  actualNumberOfThreads = runManager->GetNumberOfThreads();//This seems to not be working for bear
  if (actualNumberOfThreads<=0) actualNumberOfThreads = NumberOfThreads;
  */

  // Get MPI and Threading info
  G4int rank, size;
  if (MPI_)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get MPI rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get total number of processes (ranks)
    G4cout << "--- MPI ---\nrank: " << rank<< ", size: " <<size<< ", threads: "<<actualNumberOfThreads<<G4endl;
    // Set a unique seed for each rank
    G4Random::setTheSeed(setUniqueSeed(rank));
  }
  G4VisManager* visManager = nullptr;
  if (!MPI_)
  {
    rank = 0;
    size = 1;
    G4cout << "--- NO MPI ---\nrank: " << rank<< ", size: " <<size<< ", threads: "<<actualNumberOfThreads<<G4endl;
    visManager = new G4VisExecutive();
    visManager->Initialize();
  }


  //-- Physics list
  //The following two lines are what is used in extended/hadronic/FissionFragment
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments(true);
  G4PhysListFactory factory;
  G4String physName = "";

  G4VModularPhysicsList* physicsList = nullptr;

  char* path = std::getenv("PHYSLIST");
  if (path) { physName = G4String(path); }

  if(physName != "" && factory.IsReferencePhysList(physName))
  {
    physicsList = factory.GetReferencePhysList(physName);
    std::cout << "Got physlist "<< physName << std::endl;
  }
  else if (physName != "" && !factory.IsReferencePhysList(physName))
  {
    G4cerr << "Error: " << physName << " is not a valid physics list." << G4endl;
    return 1;
  }
  else if (physName=="")  physicsList = new PhysicsList();


  runManager->SetUserInitialization(physicsList);
  if (physName) G4cout <<"Setting up physics "<<physName<<G4endl;
  else G4cout<<"Setting up default PhysicsList"<<G4endl;

  //Initialise Additional Units
  AddUnits();
  //Set DetectorConstruction
  DetectorConstruction* fDetectorConstruction = new DetectorConstruction();
  fDetectorConstruction->SetWaterBath(noWaterBath);
  fDetectorConstruction->SetIsotope(IsotopeString);
  fDetectorConstruction->SetCasing(CasingSelection);
  fDetectorConstruction->SetAzimuthalScoring(AzimuthalScoring);
  runManager->SetUserInitialization(fDetectorConstruction);

  //-- Action Initialisation
  ActionInitialization* fActionInitialization = new ActionInitialization(rank, actualNumberOfThreads, fDetectorConstruction, FissFragments, IsotopeString, InitialNeutrons, ScoreGamma, AzimuthalScoring);
  runManager->SetUserInitialization(fActionInitialization);

  if (MPI_)
  {
    visManager = new G4VisExecutive();
    visManager->Initialize();
    g4MPISession->SessionStart();
  }
  else
  {
    G4UIExecutive* ui = nullptr;
    if (interactive)
      {
        G4cout << "Creating interactive UI session ...";
        ui = new G4UIExecutive(argc, argv);
      }
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    for (auto macro : macros)
    {
      G4String command = "/control/execute ";
      UImanager->ApplyCommand(command + macro);
    }

    if (interactive)
    {
      if (ui->IsGUI())
      {
        UImanager->ApplyCommand("/control/execute macros/ui.mac");
      }
      else
      {
        UImanager->ApplyCommand("/run/initialize");
      }
      ui->SessionStart();
      delete ui;
    }
  }

  delete visManager;
  delete g4MPI;
  delete runManager;

  G4cout << "Application successfully ended.\nBye :-)" << G4endl;

  return 0;
}

G4long setUniqueSeed(int rank)
{
  unsigned int urandomSeed;
  std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);
  if (urandom)
  {
    urandom.read(reinterpret_cast<char*>(&urandomSeed), sizeof(urandomSeed));
    urandom.close();
  }
  else
  {
    urandomSeed = 0;  // Fallback if /dev/urandom isn't available
  }

  G4long seed = time(NULL) ^ getpid() ^ (rank * 1000) ^ urandomSeed;
  return seed;
}
