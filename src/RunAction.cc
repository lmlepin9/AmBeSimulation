#include "RunAction.hh"
#include "Analysis.hh"

#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4Neutron.hh>
#include <G4AccumulableManager.hh>

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

#include <G4Threading.hh>
#include <TROOT.h>
#include <TString.h>
#include <TGraph2DErrors.h>


#include "G4AutoLock.hh"

namespace
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(int rank, int NumberOfThreads, bool FissFragments, bool NeutronTracking, bool InitialNeutrons, G4String RadioIsotope, bool ScoreGamma, bool AzimuthalScoring, bool SaveEmerging):
  G4UserRunAction(),
  fRank(rank),
  fNumberOfThreads(NumberOfThreads),
  fFissFragments(FissFragments),
  fNeutronTracking(NeutronTracking),
  fInitialNeutrons(InitialNeutrons),
  fRadioIsotope(RadioIsotope),
  fScoreGamma(ScoreGamma),
  fAzimuthalScoring(AzimuthalScoring),
  fSaveEmerging(SaveEmerging)
{
  //Messenger
  RunActionMessenger();


  fThreadid = G4Threading::G4GetThreadId();//omp_get_thread_num();
  if(fThreadid<0) fThreadid=0;
  // if (fThreadid<0 || fThreadid>fNumberOfThreads-1) G4cout << "ERROR AT THREAD NUMBER: " << fThreadid << G4endl;

  //Multithreaded protected Code
  if (true)
  {
    ROOT::EnableThreadSafety();//see https://root-forum.cern.ch/t/mutexes-when-running-inside-geant4-threads/42060
    //mutex lock for multithreaded unsafe code
        G4AutoLock l(&aMutex);
    //Hard-code diff Xsection for 9Be(a,n)12C - Part of PrimaryGenerator
    gXS_0 = new TGraph("AmBeData/GVDZ1976/GVDZ1976-intXSN0.dat");
    gXS_0->SetName("gXS_0"+fThreadid);
    gXS_0->SetTitle("gXS_0"+fThreadid);
    gXS_1 = new TGraph("AmBeData/GVDZ1976/GVDZ1976-intXSN1.dat");
    gXS_1->SetName("gXS_1"+fThreadid);
    gXS_1->SetTitle("gXS_1"+fThreadid);
    gXS_2 = new TGraph("AmBeData/GVDZ1976/GVDZ1976-intXSN2.dat");
    gXS_2->SetName("gXS_2"+fThreadid);
    gXS_2->SetTitle("gXS_2"+fThreadid);
    gXS_3 = new TGraph("AmBeData/GVDZ1976/GVDZ1976-intXSN3.dat");
    gXS_3->SetName("gXS_3"+fThreadid);
    gXS_3->SetTitle("gXS_3"+fThreadid);
    gstopping = new TGraph("AmBeData/stopping_power.txt");
    gstopping->SetName("gstopping"+fThreadid);
    gstopping->SetTitle("gstopping"+fThreadid);
    gXS_t = new TGraph("AmBeData/GVDZ1976/GVDZ1976-intXSManualSum.dat");
    gXS_t->SetName("gXS_t"+fThreadid);
    gXS_t->SetTitle("gXS_t"+fThreadid);

    /*
    //g_gs = loadAng("AmBeData/n0_ang.dat");//CM
    g_gs = loadAng("AmBeData/n0_ang_LAB.dat");
    g_gs->SetName("g_gs"+fThreadid);
    g_gs->SetTitle("g_gs"+fThreadid);
    //*/
    g_gs = loadLegendre("AmBeData/GVDZ1976/GVDZ1976-DiffN0.dat");
    if (fDBG) G4cout << "TGraphVector N0 thread" <<fThreadid <<G4endl;//DBG
    for (long unsigned int i=0; i<g_gs.size(); i++) if (g_gs[i])
    {
      g_gs[i]->SetName("g_gs"+fThreadid);
      g_gs[i]->SetTitle("g_gs"+fThreadid);
    }
    //gp = loadAng("AmBeData/n1_ang.dat");//CM
    gp = loadLegendre("AmBeData/GVDZ1976/GVDZ1976-DiffN1.dat");
    if (fDBG) G4cout << "TGraphVector N1 thread" <<fThreadid <<G4endl;//DBG
    for (long unsigned int i=0; i<gp.size(); i++) if (gp[i])
    {
      gp[i]->SetName("gp"+fThreadid);
      gp[i]->SetTitle("gp"+fThreadid);
    }
    //gpp = loadAng("AmBeData/n2_ang.dat");//CM
    gpp = loadLegendre("AmBeData/GVDZ1976/GVDZ1976-DiffN2.dat");
    if (fDBG) G4cout << "TGraphVector N2 thread " <<fThreadid <<G4endl;//DBG
    for (long unsigned int i=0; i<gpp.size(); i++) if (gpp[i])
    {
      gpp[i]->SetName("gpp"+fThreadid);
      gpp[i]->SetTitle("gpp"+fThreadid);
    }
    gppp = loadLegendre("AmBeData/GVDZ1976/GVDZ1976-DiffN3.dat");
    if (fDBG) G4cout << "TGraphVector N3 thread " <<fThreadid <<G4endl;//DBG
    for (long unsigned int i=0; i<gppp.size(); i++) if (gppp[i])
    {
      gppp[i]->SetName("gppp"+fThreadid);
      gppp[i]->SetTitle("gppp"+fThreadid);
    }

  }


  //Analysis via Root TTrees
  if (fRadioIsotope=="USF") FileOutNameRank = fRadioIsotope + FileOutName+"-N"+fRank;
  else if (fRadioIsotope.length()>=6 and fRadioIsotope.substr(0,6)=="Single") FileOutNameRank = fRadioIsotope + FileOutName+"-N"+fRank;
  else FileOutNameRank = fRadioIsotope.substr(3,5) + "Be" + FileOutName+"-N"+fRank;
  fout.resize(fNumberOfThreads);
  tree.resize(fNumberOfThreads);
  if(fout[fThreadid]==0)
  {
    G4cout<<"THREAD: "<<FileOutNameRank+"-T"+fThreadid+".root"<<G4endl;
    fout[fThreadid]=new TFile(FileOutNameRank+"-T"+fThreadid+".root","RECREATE");
    tree[fThreadid] = new TTree(TString("FullBath"), "FullBath");//TTree name
    tree[fThreadid]->SetDirectory(fout[fThreadid]);


    if (true)
    {
      tree[fThreadid]->Branch("EmergingNeutrons",&fEmergingNeutrons);

      //tree[fThreadid]->Branch("secondaryElec",&fbranchSecondaryElectron);

      if (fScoreGamma)//Score secondary gammas to compare with hpGe
      {
        tree[fThreadid]->Branch("secondaryGamma",&fbranchSecondaryGamma);
        tree[fThreadid]->Branch("secondaryGammaEmerging",&fbranchSecondaryGammaEmerging);
      }
    }
    if (fInitialNeutrons)
    {
      tree[fThreadid]->Branch("initN",&fbranchEInitN);
    }
    if (fNeutronTracking)
    {
      tree[fThreadid]->Branch("secondaryNeut",&fbranchSecondaryNeutron);
      tree[fThreadid]->Branch("secondaryNeutEmerging",&fbranchSecondaryNeutronEmerging);
    }
    if (fFissFragments)
    {
      tree[fThreadid]->Branch("FissProds",&fbranchFissionProduct);
      tree[fThreadid]->Branch("FissProdMass",&fbranchFissionProductMass);
      tree[fThreadid]->Branch("FissIon",&fbranchFissIon);
      tree[fThreadid]->Branch("FissNeut",&fbranchFissNeut);
      tree[fThreadid]->Branch("FissNeutEmerging",&fbranchFissNeutEmerging);
    }
    if (fAzimuthalScoring)
    {
      tree[fThreadid]->Branch("fbranchNEmissionSpec",&fbranchNEmissionSpec);
      tree[fThreadid]->Branch("fbranchNEmissionSpecVer",&fbranchNEmissionSpecVer);
    }
  }


  // We borrow this structure to store our emerging particles 

  if(fSaveEmerging){
      if (fRadioIsotope=="USF") EmergingFileOutNameRank = fRadioIsotope + EmergingFileOutName+"-N"+fRank;
      else if (fRadioIsotope.length()>=6 and fRadioIsotope.substr(0,6)=="Single") EmergingFileOutNameRank = fRadioIsotope + EmergingFileOutName+"-N"+fRank;
      else EmergingFileOutNameRank = fRadioIsotope.substr(3,5) + "Be" + EmergingFileOutName+"-N"+fRank;
      fEmergingOut.resize(fNumberOfThreads);
      fEmergingTree.resize(fNumberOfThreads);

      if(fEmergingOut[fThreadid]==0)
      {
        G4cout<<"THREAD: "<<EmergingFileOutNameRank+"-T"+fThreadid+".root"<<G4endl;
        fEmergingOut[fThreadid]=new TFile(EmergingFileOutNameRank+"-T"+fThreadid+".root","RECREATE");
        fEmergingTree[fThreadid] = new TTree(TString("EmergingParticles"), "EmergingParticles");//TTree name
        fEmergingTree[fThreadid]->SetDirectory(fEmergingOut[fThreadid]);

        fEmergingTree[fThreadid]->Branch("TrackId",&fbranchEmergingId);
        fEmergingTree[fThreadid]->Branch("ParentId",&fbranchEmergingParentId);
        fEmergingTree[fThreadid]->Branch("PDG",&fbranchEmergingPDG);
        fEmergingTree[fThreadid]->Branch("Vertex",&fbranchEmergingPos);
        fEmergingTree[fThreadid]->Branch("Momentum", &fbranchEmergingP);
        fEmergingTree[fThreadid]->Branch("Process", &fbranchEmergingProcess);
        
      }
  }





}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  G4cout<<"Trying to close thread "<<fThreadid<<" for "<<fout[fThreadid]<<G4endl;
  if (fout[fThreadid]!=0) fout[fThreadid]->Close();
  if (fEmergingOut[fThreadid]!=0) fEmergingOut[fThreadid]->Close();

  if (gXS_0) delete gXS_0;
  if (gXS_1) delete gXS_1;
  if (gXS_2) delete gXS_2;
  if (gXS_3) delete gXS_3;
  if (gstopping) delete gstopping;
  if (gXS_t) delete gXS_t;
  //if (g_gs) delete g_gs;
  //if (gp) delete gp;
  //if (gpp) delete gpp;
  //if (gppp) delete gppp;

  if (IsMaster())//merge into one file per rank
  {
    //Merge root files
    G4String hadd_cmd = "hadd -f " + FileOutNameRank + ".root " + FileOutNameRank + "-T*.root";
    G4cout << "Running hadd command" << hadd_cmd.c_str() << G4endl;
    int haddFiles = system((hadd_cmd).c_str());

    if (haddFiles==0)
    {
      int removeFiles = system(("rm -vf "+FileOutNameRank+"-T*.root").c_str());
      if (removeFiles==0) G4cout << "Files merged and removed successfully" << G4endl;
      else G4cout << "Files merged but NOT removed" << G4endl;
    }
    else G4cout << "Files NOT merged - check manually" << G4endl;
  }


  // Merge emerging particles files 
  if(IsMaster() && fSaveEmerging){
    //Merge root files
    G4String hadd_cmd = "hadd -f " + EmergingFileOutNameRank + ".root " + EmergingFileOutNameRank + "-T*.root";
    G4cout << "Running hadd command" << hadd_cmd.c_str() << G4endl;
    int haddFiles = system((hadd_cmd).c_str());

    if (haddFiles==0)
    {
      int removeFiles = system(("rm -vf "+EmergingFileOutNameRank+"-T*.root").c_str());
      if (removeFiles==0) G4cout << "Files merged and removed successfully" << G4endl;
      else G4cout << "Files merged but NOT removed" << G4endl;
    }
    else G4cout << "Files NOT merged - check manually" << G4endl;

  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;


  G4cout<<"Finished thread: "<<fThreadid<<G4endl;
  //G4cout<<"Run ID: "<<thr<<G4endl;

  ///*
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events " << G4endl;
  }

  //DO NOT replace with TreeWriteMidRun
  //Close File
  if(fout[fThreadid]!=0 && fThreadid!=-1)
  {
    //G4cout<<"FILL "<<fThreadid<<"\t"<<tree[fThreadid]<<"\t"<<fout[fThreadid]<<G4endl;
    //tree[fThreadid]->Fill();//Don't fill here as it would add an extra event
    fout[fThreadid]->cd();
    tree[fThreadid]->Print();
    G4cout<<"Write"<<G4endl;
    tree[fThreadid]->Write();
    G4cout<<"CLOSE"<<G4endl;
  }

  if(fEmergingOut[fThreadid]!=0 && fThreadid!=-1 && fSaveEmerging){
    fEmergingOut[fThreadid]->cd();
    fEmergingTree[fThreadid]->Print();
    G4cout<<"Write"<<G4endl;
    fEmergingTree[fThreadid]->Write();
    G4cout<<"CLOSE"<<G4endl;
  }


  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TreeFill(G4int eventNumber)//end of event
{
  fEventNo = eventNumber;//this increases also for empty events, which are not recorded
  if (fout[fThreadid]!=0)
  {
    tree[fThreadid]->Fill();
    //G4cout<<"FILL "<<fThreadid<<"\t"<<tree[fThreadid]<<"\t"<<fout[fThreadid]<<G4endl;
  }

  if(fEmergingOut[fThreadid]!=0 && fSaveEmerging){
    fEmergingTree[fThreadid]->Fill();
  }
  //reset variables
  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ClearBranches()
{
  //non-vecotr data can be simply assigned and overwritten
  fbranchEInitN = -1;
  fEventNo = -1;
  fbranchValid = false;

  fbranchSecondaryNeutron.clear();
  fbranchSecondaryNeutronEmerging.clear();
  fbranchSecondaryElectron.clear();
  fbranchSecondaryElectronEmerging.clear();
  fbranchSecondaryGamma.clear();
  fbranchSecondaryGammaEmerging.clear();

  fbranchFissionProduct.clear();
  fbranchFissionProductMass.clear();

  fEmergingNeutrons.clear();

  fbranchFissIon.clear();
  fbranchFissNeut.clear();
  fbranchFissNeutEmerging.clear();

  fbranchNEmissionSpec.clear();
  fbranchNEmissionSpecVer.clear();

  fbranchEmergingId.clear();
  fbranchEmergingParentId.clear();
  fbranchEmergingPDG.clear();
  fbranchEmergingPos.clear();
  fbranchEmergingP.clear(); 
  fbranchEmergingProcess.clear(); 


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RecordSecondaries(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma)
{
  fbranchSecondaryNeutron = neutron;
  fbranchSecondaryElectron = electron;
  fbranchSecondaryGamma = gamma;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RecordSecondariesEmerging(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma)
{
  fbranchSecondaryNeutronEmerging = neutron;
  fbranchSecondaryElectronEmerging = electron;
  fbranchSecondaryGammaEmerging = gamma;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::loadFile(TString filename, std::vector<double>& Z_array, std::vector<double>& X_array, std::vector<double>& Y_array, std::vector<double>& ex_array, std::vector<double>& ey_array)
{
  std::ifstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::string line;
  while (std::getline(file, line))
  {
    std::stringstream ss(line);
    double Z, X, Y, ex, ey;

    if (ss >> Z >> X >> Y >> ex >> ey)
    {
      Z_array.push_back(Z);
      X_array.push_back(X);
      Y_array.push_back(Y);
      ex_array.push_back(ex);
      ey_array.push_back(ey);
    }
    else
    {
      std::cerr << "Warning: Could not parse line: " << line << std::endl;
    }
  }

  file.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RunActionMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/run/", "Run control");
  fMessenger->DeclareProperty("Debug",fDBG);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TGraph2DErrors* RunAction::loadAng(TString filename)
{
    std::vector<double> Z_array, X_array, Y_array, ex_array, ey_array, ez_array;
    loadFile(filename, Z_array, X_array, Y_array, ex_array, ey_array);
    for (long unsigned int i=0; i<Z_array.size(); i++) ez_array.push_back(0);

    // Display loaded data
    TGraph2DErrors *tg = new TGraph2DErrors(Z_array.size(),&Z_array[0],&X_array[0],&Y_array[0],&ez_array[0],&ex_array[0],&ey_array[0]);
    return tg;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<TGraph*> RunAction::loadLegendre(TString filename)
{
  std::vector<Double_t> n_arr, X_arr, Y_arr;

  //Open and read file
  std::ifstream file(filename);
  if (!file.is_open())
  {
    G4cerr << "Error: Could not open file " << filename << G4endl;
    //return;
  }

  G4cout << "Opening legendre pol file" <<G4endl;

  std::string line;
  while (std::getline(file, line))
  {
    std::stringstream ss(line);
    double n, X, Y;

    if (ss >> n >> X >> Y)
    {
      n_arr.push_back(n);
      X_arr.push_back(X);
      Y_arr.push_back(Y);
    }
    else
    {
      G4cerr << "Warning: Could not parse line: " << line << G4endl;
    }
  }
  file.close();

  if (fDBG) G4cout << "RunAction::loadLegendre - file read successfully"<<G4endl;

  Int_t max_n = static_cast<Int_t>(*std::max_element(n_arr.begin(), n_arr.end()));

  //create TGraphs
  std::vector<TGraph*> tg;
  tg.reserve(max_n + 1);

  for (Int_t i=0; i<=max_n; i++)
  {
    std::vector<Double_t> x, y;
    for (std::size_t j = 0; j < n_arr.size(); ++j)
    {
      if (static_cast<Int_t>(n_arr[j]) == i)
      {
        x.push_back(X_arr[j]);
        y.push_back(Y_arr[j]);
      }
    }

    if (!x.empty())
    {
      tg.push_back(new TGraph(x.size(), x.data(), y.data()));

      // --- Debug save as PNG ---
      //TCanvas c(TString::Format("c_n%d", i), TString::Format("c_n%d", i), 800, 600);
      //tg[i]->Draw("ALP");
      //c.SaveAs(Form("tg_n%d.png", i));
    }
    else
    {
      tg[i] = nullptr; // optional, keeps vector aligned
    }
  }
  return tg;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
