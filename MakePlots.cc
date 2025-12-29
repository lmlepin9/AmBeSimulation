/*
 * Root macro?
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <vector>
#include <array>
#include <math.h>
//Unix folder creation
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//ROOT libraries
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TSystem.h>

#include <ROOT/RDataFrame.hxx>

#include <Math/Vector3D.h>
//#include <Math/Vector3Dfwd.h>
//From Geant4
#include <CLHEP/Vector/ThreeVector.h>


//constants
Double_t RunTime = 100;//s

int MakePlots(TFile* _file0);
int MakePlots(TString inname="");

Bool_t fEmergingGammas=true;
Bool_t fFission=true;
Bool_t isotropic=false;

//code
int MakePlots(TFile* _file0)
{
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel

  ROOT::RDataFrame tree_FullBath("FullBath", _file0); // Interface to TTree and TChain

  std::vector<TString> BranchList = {"FissIon", "secondaryGamma", "secondaryGammaEmerging", "EmergingNeutrons", "PerpNSpectrumEmerging", "VertNSpectrumEmerging", "secondaryNeut", "secondaryNeutEmerging", "initN", "FissNeut", "FissNeutEmerging", };

  if (!fEmergingGammas)
  {
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "secondaryGamma"), BranchList.end());
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "secondaryGammaEmerging"), BranchList.end());
  }
  if (!fFission)
  {
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "FissIon"), BranchList.end());
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "FissNeut"), BranchList.end());
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "FissNeutEmerging"), BranchList.end());
  }
  if (!isotropic)
  {
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "PerpNSpectrumEmerging"), BranchList.end());
    BranchList.erase(std::remove(BranchList.begin(), BranchList.end(), "VertNSpectrumEmerging"), BranchList.end());
  }

  TString filename = _file0->GetName();
  TString date = filename(0,4)+filename(5,2)+filename(8,2);
  TString filename_base = filename(11,filename.Length()-5-11);
  filename_base.ReplaceAll("-","_");
  TString filename_new = filename_base+"_"+date+"_";

  TCanvas c("c","x hist");
  ROOT::RDF::RResultPtr<TH1D> myHisto;

  unsigned int i = 0;
  if (fFission)
  {
    myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 150, 19., 199.}, BranchList[i]);
    myHisto->Draw();
    c.SaveAs(filename_new+BranchList[i]+".C");
    i++;
  }
  if (fEmergingGammas)
  {
    myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 1000, 0., 10.}, BranchList[i]);
    myHisto->Draw();
    c.SaveAs(filename_new+BranchList[i]+".C");
    i++;
    myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 1000, 0., 10.}, BranchList[i]);
    myHisto->Draw();
    c.SaveAs(filename_new+BranchList[i]+".C");
    i++;
  }

  for (; i<BranchList.size(); i++)
  {
    //myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 120, 0., 24.}, BranchList[i]);
    myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 480, 0., 24.}, BranchList[i]);
    //myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 1000, 0., 10.}, BranchList[i]);
    myHisto->Draw();
    c.SaveAs(filename_new+BranchList[i]+".C");
  }
  //*/



  return 0;
}

int MakePlots(TString inname)
{
 ///---OPEN AND READ FILE---///
  TString filename;
  if (inname=="")
  {
    std::cout << "Enter filename to analyse" << "\n";
    std::cin >> filename;
  }
  else
  {
    filename=inname;
  }
  std::cout <<"Now analysing: " << filename << "\n";//MSG
  //add root extension if missing
  if (filename(filename.Length()-5,filename.Length()) != ".root") filename=filename+".root";

  //initialise TTree read
  if(gSystem->AccessPathName(filename))
  {
    std::cout << "File does not exist" << std::endl;//MSG
    return 1;
  }
  TFile* file = new TFile(filename, "read");

  return MakePlots(file);
}
