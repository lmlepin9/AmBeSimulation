#include "DetectorConstruction.hh"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4MultiUnion.hh>

#include <G4SDManager.hh>
#include <G4SDParticleFilter.hh>

#include <G4MultiFunctionalDetector.hh>
#include <G4PSEnergyDeposit.hh>
#include <G4PSCellFlux.hh>
#include <G4PSDoseDeposit.hh>

#include <sstream>

//using namespace std;


/*
DetectorConstruction::DetectorConstruction(G4int worldType, G4String IsotopeString)
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fworldType(worldType),
  fRadioIsotope(IsotopeString)
//*/
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  fDebug = false;

  fCasingSelection = 0;
  DetectorMessenger();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4int ncomponents, natoms;
  G4bool checkOverlaps = true;
  G4NistManager* nist = G4NistManager::Instance();

  //--Elements
  G4Element* C = nist->FindOrBuildElement("C");
  G4Element* O = nist->FindOrBuildElement("O");
  G4Element* Si = nist->FindOrBuildElement("Si");
  G4Element* P = nist->FindOrBuildElement("P");
  G4Element* S = nist->FindOrBuildElement("S");
  G4Element* Cr = nist->FindOrBuildElement("Cr");
  G4Element* Mn = nist->FindOrBuildElement("Mn");
  G4Element* Fe = nist->FindOrBuildElement("Fe");
  G4Element* Ni = nist->FindOrBuildElement("Ni");
  G4Element* Mo = nist->FindOrBuildElement("Mo");


  //--Materials
  G4Material* matBe = nist->FindOrBuildMaterial("G4_Be");
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* matVacuum = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* matWater = nist->FindOrBuildMaterial("G4_WATER");
  G4cout<<"WATER TEMP: "<<matWater->GetTemperature()<<G4endl;

  G4Material* Steel304 = new G4Material("Steel304", 8*g/cm3, ncomponents=3);
  Steel304->AddElement(Fe, 72*perCent);
  Steel304->AddElement(Cr, 18*perCent);
  Steel304->AddElement(Ni,  10*perCent);

  //Steel: 02X17H14M2
  //from https://www.sciencedirect.com/science/article/pii/S0022311504002521?via=ihub
  G4Material* matInox = new G4Material("Stainless-Steel", 7.917*g/cm3, ncomponents=9);
  matInox->AddElement(Fe, 64.653*perCent);
  matInox->AddElement(C, 0.02*perCent);
  matInox->AddElement(Si, 0.4*perCent);
  matInox->AddElement(Mn, 1.5*perCent);
  matInox->AddElement(P, 0.017*perCent);
  matInox->AddElement(S, 0.01*perCent);
  matInox->AddElement(Cr, 17.0*perCent);
  matInox->AddElement(Ni, 14.0*perCent);
  matInox->AddElement(Mo, 2.4*perCent);
  //*/

  //D2O
  G4Isotope* isoDeut = new G4Isotope("isoDeut", 1, 2, 2. * g/mole);
  G4Element* elDeut = new G4Element("elDeut", "elDeut", 1);
  elDeut->AddIsotope(isoDeut, 100.*perCent);

  G4Material* matD2O = new G4Material("D2O", 1.1056 *g/cm3, ncomponents=2);
  matD2O->AddElement(elDeut, 2);
  matD2O->AddElement(O, 1);


  //Am241
  G4Isotope* isoAm = new G4Isotope("isoAm", 95, 241, 241. * g/mole);
  G4Element* elAm = new G4Element("Am241", "Am241", 1);
  elAm->AddIsotope(isoAm, 100.*perCent);

  //Am243
  G4Isotope* isoAm243 = new G4Isotope("isoAm243", 95, 243, 243. * g/mole);
  G4Element* elAm243 = new G4Element("Am243", "Am243", 1);
  elAm243->AddIsotope(isoAm243, 100.*perCent);

  //Pu239
  G4Isotope* isoPu = new G4Isotope("isoAm", 94, 239, 239. * g/mole);
  G4Element* elPu = new G4Element("Pu239", "Pu239", 1);
  elPu->AddIsotope(isoPu, 100.*perCent);

  //Am2O3
  G4Material* matAm2O3 = new G4Material("Am2O3", 11.77*g/cm3, ncomponents=2);
  matAm2O3->AddElement(elAm,natoms=2);
  matAm2O3->AddElement(O,natoms=3);

  //AmO2
  G4Material* matAmO2 = new G4Material("AmO2", 11.68*g/cm3, ncomponents=2);
  matAmO2->AddElement(elAm,natoms=1);
  matAmO2->AddElement(O,natoms=2);

  //PuO2
  G4Material* matPuO2 = new G4Material("PuO2", 11.5*g/cm3, ncomponents=2);
  matPuO2->AddElement(elPu,natoms=1);
  matPuO2->AddElement(O,natoms=2);

  //AmBe
  G4double PerBe = .890326;
  G4double PerAmO2 = 1-PerBe;
  G4double densityAmBe = (1.848*PerBe+11.68*PerAmO2);//2.9265704g/cm3
  G4double densityPuBe = (1.848*PerBe+11.5*PerAmO2);//2.9265704g/cm3

  G4Material* matAmBe = new G4Material("AmBe", densityAmBe*g/cm3, ncomponents=2);
  matAmBe->AddMaterial(matBe, PerBe*100*perCent);
  if (fRadioIsotope=="241Am")
    matAmBe->AddMaterial(matAmO2, PerAmO2*100*perCent);
  else if (fRadioIsotope=="239Pu")
    matAmBe->AddMaterial(matPuO2, PerAmO2*100*perCent);
  else matAmBe->AddMaterial(matAmO2, PerAmO2*100*perCent);


  //matAmBe->AddMaterial(matVacuum, PerAmO2*100*perCent);//TEST ONLY


  //--Colours
  //red
  G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());
  red->SetVisibility(true);
  //red->SetForceSolid(true);
  //yellow
  G4VisAttributes* yellow = new G4VisAttributes(G4Colour::Yellow());
  yellow->SetVisibility(true);
  //yellow->SetForceSolid(true);
  //blue
  G4VisAttributes* blue = new G4VisAttributes(G4Colour::Blue());
  blue->SetVisibility(true);
  //blue->SetForceSolid(true);
  //green
  G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());
  green->SetVisibility(true);
  green->SetForceSolid(true);
  //magenta
  G4VisAttributes* magenta = new G4VisAttributes(G4Colour::Magenta());
  magenta->SetVisibility(true);
  magenta->SetForceSolid(true);
  //cyan
  G4VisAttributes* cyan = new G4VisAttributes(G4Colour::Cyan());
  cyan->SetVisibility(true);
  cyan->SetForceSolid(true);


  //--Measurements
  //absorber
  //G4double AbsorberRadius = 17.6*mm/2.;
  G4double AbsorberRadius = 17.5*mm/2.;//X3 Specs
  G4double AbsorberLength = 17.5*mm;
  G4double ContainThickness = 2.4*mm;

  /*
  //absorber -- RUS/6137/S-96 IBN-241-15-4
  G4double AbsorberRadius = 17.4*mm/2.;
  G4double AbsorberLength = 17.6*mm;
  G4double ContainThicknessR = 2.54*mm;
  G4double ContainThicknessL = 2.4*mm;
  //*/


  //--World
  //G4double WorldSizeXY = 50.*cm *3.*10;
  //G4double WorldSizeZ  = 51.*cm *3.*10;
  G4double WorldSizeXY = 50*cm*3.*10;
  G4double WorldSizeZ  = 50*cm*3.*10;

  if (fworldType==-1)
  {
    WorldSizeZ = 500 * cm;
    WorldSizeXY = 500 * cm;
  }


  // compute dimensions
  G4double ContainRadius = AbsorberRadius + ContainThickness;
  G4double ContainLength = AbsorberLength + 2*ContainThickness;



  //-- GEOMETRY
  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(true);

  //--World
  G4Tubs* WorldSolid = new G4Tubs("WorldSolid", 0., WorldSizeXY, WorldSizeZ, 0., 360.*deg);
  G4LogicalVolume* WorldLog = nullptr;
  //WorldLog= new G4LogicalVolume(WorldSolid, matAir, "WorldLogical");
  if (fworldType==2)
    WorldLog = new G4LogicalVolume(WorldSolid, matWater, "WorldLogical");
  else if (fworldType==3)
    WorldLog = new G4LogicalVolume(WorldSolid, matD2O, "WorldLogical");
  else
    WorldLog = new G4LogicalVolume(WorldSolid, matAir, "WorldLogical");
  WorldLog->SetVisAttributes(visAttr);
  G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), WorldLog, "World", 0, false, 0);


  //X3 Casing
  G4RotationMatrix rotm = G4RotationMatrix();
  G4RotationMatrix rotm180 = G4RotationMatrix();
  rotm180.rotateX(180*deg);

  //Item3 -- Addition
  G4double IT3_r1 = 19.8/2 *mm;
  G4double IT3_r2 = 17.5/2 *mm;
  G4double IT3_r3 = 19.3/2 *mm;
  G4double IT3_h1 = 21.7 *mm;
  G4double IT3_h2 = 20.5 *mm;
  G4double IT3_h3 = 2.0 *mm;

  G4VSolid* IT3BaseSolid = new G4Tubs("IT3BaseSolid", 0., IT3_r1, (IT3_h1-IT3_h2)/2, 0, 360*deg);
  //G4VSolid* IT3CoreSolid = new G4Tubs("IT3CoreSolid", IT3_r2, IT3_r1, (IT3_h2-IT3_h3)/2, 0, 360*deg);
  //Use 0 as inner radius as AmBe needs to be a daughter volume of the X3 casing
  G4VSolid* IT3CoreSolid = new G4Tubs("IT3CoreSolid", 0, IT3_r1, (IT3_h2-3.0*mm)/2, 0, 360*deg);
  G4VSolid* IT3CoreTSolid = new G4Tubs("IT3CoreTSolid", IT3_r2, IT3_r1, (1.0*mm)/2, 0, 360*deg);
  G4VSolid* IT3TopSolid = new G4Tubs("IT3TopSolid", IT3_r2, IT3_r3, (IT3_h3)/2, 0, 360*deg);

  G4MultiUnion* IT3Solid = new G4MultiUnion("IT3Solid");
  IT3Solid->AddNode(*IT3BaseSolid, G4Transform3D(rotm,G4ThreeVector(0.,0.,-IT3_h1/2 + (IT3_h1-IT3_h2)/2)));
  IT3Solid->AddNode(*IT3CoreSolid, G4Transform3D(rotm,G4ThreeVector(0.,0.,-IT3_h1/2 + (IT3_h1-IT3_h2) + (IT3_h2-3.0*mm)/2)));
  IT3Solid->AddNode(*IT3CoreTSolid, G4Transform3D(rotm,G4ThreeVector(0.,0.,-IT3_h1/2 + (IT3_h1-IT3_h2) + (IT3_h2-3.0*mm) + (1.0*mm)/2)));
  IT3Solid->AddNode(*IT3TopSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., IT3_h1/2 - IT3_h3/2)));
  IT3Solid->Voxelize();

  G4LogicalVolume* IT3Log = new G4LogicalVolume(IT3Solid, Steel304, "IT3Logical");
  //IT3Log->SetVisAttributes(magenta);
  //G4PVPlacement* IT3Phys = new G4PVPlacement(0, G4ThreeVector(), IT3Log, "IT3", WorldLog, false, 0, checkOverlaps);


  //Item4 -- Addition
  G4double IT4_r1 = 15.6/2 *mm;
  G4double IT4_h1 = 3.0 *mm;
  G4double IT4_h2 = 1.0 *mm;

  G4VSolid* IT4BottomSolid = new G4Tubs("IT4BottomSolid", 0, IT3_r2, (IT4_h1-IT4_h2)/2, 0, 360*deg);
  G4VSolid* IT4TopSolid = new G4Tubs("IT4TopSolid", IT4_r1, IT3_r2, (IT4_h2)/2, 0, 360*deg);

  G4MultiUnion* IT4Solid = new G4MultiUnion("IT4Solid");
  // Bottom: z from -IT4_h1/2 to -IT4_h1/2+(IT4_h1-IT4_h2)
  IT4Solid->AddNode(*IT4BottomSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -(IT4_h1-IT4_h2)/2)));
  // Top: z from -IT4_h1/2+(IT4_h1-IT4_h2) to IT4_h1/2
  IT4Solid->AddNode(*IT4TopSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., + IT4_h2/2)));
  IT4Solid->Voxelize();
  G4LogicalVolume* IT4Log = new G4LogicalVolume(IT4Solid, Steel304, "IT4Logical");
  //IT4Log->SetVisAttributes(cyan);
  //G4PVPlacement* IT4Phys = new G4PVPlacement(0, G4ThreeVector({0.,0.,100.*mm}), IT4Log, "IT4", WorldLog, false, 0, checkOverlaps);


  //Item2
  G4double IT2_r1 = 20.0/2 *mm;
  G4double IT2_r2 = 17.3/2 *mm;
  G4double IT2_r3 = 10.0/2 *mm;
  G4double IT2_r4 = 18.2/2 *mm;
  G4double IT2_rM6 = 6.0/2 *mm;
  G4double IT2_h1 = 6.0 *mm;
  G4double IT2_h2 = 4.0 *mm;
  G4double IT2_h3 = 2.0 *mm;
  G4double IT2_h4 = 0.5 *mm;

  G4VSolid* IT2BottomSolid = new G4Tubs("IT2BottomSolid", 0, IT2_r1, IT2_h3/2, 0, 360*deg);
  G4VSolid* IT2ExtBottomSolid = new G4Tubs("IT2ExtBottomSolid", IT2_r2, IT2_r1, IT2_h3/2, 0, 360*deg);
  G4VSolid* IT2ExtTopSolid = new G4Tubs("IT2ExtTopSolid", IT2_r4, IT2_r1, IT2_h3/2, 0, 360*deg);
  G4VSolid* IT2IntSolid = new G4Tubs("IT2IntSolid", IT2_rM6, IT2_r3, (IT2_h1-IT2_h3-IT2_h4)/2, 0, 360*deg);

  G4MultiUnion* IT2Solid = new G4MultiUnion("IT2Solid");
  // Bottom base: z from -IT2_h1/2 to -IT2_h1/2+IT2_h3
  IT2Solid->AddNode(*IT2BottomSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT2_h1/2 + IT2_h3/2)));
  // ExtBottom ring: z from -IT2_h1/2+IT2_h3 to -IT2_h1/2+2*IT2_h3
  IT2Solid->AddNode(*IT2ExtBottomSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT2_h1/2 + IT2_h3 + IT2_h3/2)));
  // ExtTop ring: z from -IT2_h1/2+2*IT2_h3 to -IT2_h1/2+3*IT2_h3
  IT2Solid->AddNode(*IT2ExtTopSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT2_h1/2 + 2*IT2_h3 + IT2_h3/2)));
  // Internal cylinder: z from -IT2_h1/2+IT2_h3 to IT2_h1/2-IT2_h2-IT2_h4
  IT2Solid->AddNode(*IT2IntSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT2_h1/2 + IT2_h3 + (IT2_h1-IT2_h3-IT2_h4)/2)));
  IT2Solid->Voxelize();
  //G4LogicalVolume* IT2Log = new G4LogicalVolume(IT2Solid, Steel304, "IT2Logical");
  //IT2Log->SetVisAttributes(green);
  //G4PVPlacement* IT2Phys = new G4PVPlacement(0, G4ThreeVector({0.,0.,-100.*mm}), IT2Log, "IT2", WorldLog, false, 0, checkOverlaps);


  //Item1
  G4double IT1_r1 = 22.4*mm/2;
  G4double IT1_r2 = 21.8*mm/2;
  G4double IT1_r3 = 20.0*mm/2;
  G4double IT1_rSlot = 2.0*mm/2;
  G4double IT1_h1 = 31.0*mm;
  G4double IT1_h2 = 28.0*mm;
  G4double IT1_h3 = 1.0*mm;
  G4double IT1_h4 = 2.0*mm;

  G4VSolid* IT1BaseBottomSolid = new G4Tubs("IT1BaseBottomSolid", IT1_rSlot, IT1_r1, (IT1_h3)/2, 0, 360*deg);
  G4VSolid* IT1BaseTopSolid = new G4Tubs("IT1BaseTopSolid", 0, IT1_r1, (IT1_h1-IT1_h2-IT1_h3)/2, 0, 360*deg);
  G4VSolid* IT1CylBottomSolid = new G4Tubs("IT1CylBottomSolid", 0., IT1_r1, (IT1_h2-IT2_h1)/2, 0, 360*deg);
  G4VSolid* IT1CylBottomTSolid = new G4Tubs("IT1CylBottomTSolid", IT1_r3, IT1_r1, (4.0*mm)/2, 0, 360*deg);
  G4VSolid* IT1CylTopSolid = new G4Tubs("IT1CylTopSolid", IT1_r3, IT1_r2, (IT1_h4)/2, 0, 360*deg);

  G4MultiUnion* IT1Solid = new G4MultiUnion("IT1Solid");
  // BaseBottom with slot: z from -IT1_h1/2 to -IT1_h1/2+IT1_h3
  IT1Solid->AddNode(*IT1BaseBottomSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT1_h1/2 + IT1_h3/2)));
  // BaseTop solid disc: z from -IT1_h1/2+IT1_h3 to -IT1_h1/2+IT1_h1-IT1_h2
  IT1Solid->AddNode(*IT1BaseTopSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT1_h1/2 + IT1_h3 + (IT1_h1-IT1_h2-IT1_h3)/2)));
  // CylBottom hollow: z from -IT1_h1/2+IT1_h1-IT1_h2 to IT1_h1/2-IT1_h4
  IT1Solid->AddNode(*IT1CylBottomSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT1_h1/2 + (IT1_h1-IT1_h2) + (IT1_h2-IT2_h1)/2)));
  IT1Solid->AddNode(*IT1CylBottomTSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., -IT1_h1/2 + (IT1_h1-IT1_h2) + (IT1_h2-IT2_h1) + (4.0*mm)/2)));
  // CylTop narrow hollow: z from IT1_h1/2-IT1_h4 to IT1_h1/2
  IT1Solid->AddNode(*IT1CylTopSolid, G4Transform3D(rotm,G4ThreeVector(0.,0., IT1_h1/2 - IT1_h4/2)));
  IT1Solid->Voxelize();
  //G4LogicalVolume* IT1Log = new G4LogicalVolume(IT1Solid, Steel304, "IT1Logical");
  //IT1Log->SetVisAttributes(green);
  //G4PVPlacement* IT1Phys = new G4PVPlacement(0, G4ThreeVector({0.,100.*mm,0.}), IT1Log, "IT1", WorldLog, false, 0, checkOverlaps);

  //X3 Casing
  G4MultiUnion* X3Solid = new G4MultiUnion("X3Solid");
  G4double X3ShiftZ = +0.9*mm;//it moves the X3, not the AmBe
  //X3Solid->AddNode(*IT3Solid, G4Transform3D(rotm, G4ThreeVector(0, 0, X3ShiftZ)));
  //X3Solid->AddNode(*IT4Solid, G4Transform3D(rotm, G4ThreeVector(0, 0, X3ShiftZ+(IT3_h1-3.0)*mm/2)));
  //X3Solid->AddNode(*IT2Solid, G4Transform3D(rotm180, G4ThreeVector(0, 0, X3ShiftZ+(IT3_h1+IT2_h1)*mm/2)));//BAD ORIENTATION
  X3Solid->AddNode(*IT2Solid, G4Transform3D(rotm, G4ThreeVector(0, 0, X3ShiftZ+(IT3_h1+IT2_h1)*mm/2)));
  X3Solid->AddNode(*IT1Solid, G4Transform3D(rotm, G4ThreeVector(0, 0, X3ShiftZ+(IT3_h1-31.0)*mm/2+6.0*mm)));
  X3Solid->Voxelize();
  //ContainLog->SetVisAttributes(yellow);
  //G4PVPlacement* X3Phys = new G4PVPlacement(0, G4ThreeVector({0.,0.*mm,0.}), ContainLog, "Container", WorldLog, false, 0, checkOverlaps);

  //AmBe Container selection
  // Logical for an inner vacuum when using the X3 casing
  G4LogicalVolume* X3VacLog = nullptr;

  if (fCasingSelection==0)
  {
    ContainSolid = new G4Tubs("ContainerSolid", 0., ContainRadius, 0.5*ContainLength, 0., 360.*deg);
    ContainLog = new G4LogicalVolume(ContainSolid, matInox, "ContainerLogical");
  }
  else if (fCasingSelection==1)
  {
    ContainLog = new G4LogicalVolume(X3Solid, Steel304, "ContainerLogical");
    //ContainLog = new G4LogicalVolume(X3Solid, matInox, "ContainerLogical");

    G4VSolid* X3Vacuum = new G4Tubs("X3VacuumSolid", 0., IT1_r1-0.5*mm, 0.5*(IT1_h2 - IT2_h1), 0., 360.*deg);
    X3VacLog = new G4LogicalVolume(X3Vacuum, matVacuum, "X3VacuumLogical");

    // Place vacuum at center of X3Solid coordinate system (between IT1 and IT2)
    G4PVPlacement* X3VacPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, X3ShiftZ), X3VacLog, "X3Vacuum", ContainLog, false, 0, checkOverlaps);

    // Place IT3, IT4 inside vacuum (positions relative to vacuum center)
    // IT3 centered in vacuum, IT4 sits on top of IT3
    G4PVPlacement* IT3Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), IT3Log, "IT3", X3VacLog, false, 0, checkOverlaps);
    G4PVPlacement* IT4Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, IT3_h1/2 + IT4_h1/2 - 2.5*mm), IT4Log, "IT4", X3VacLog, false, 0, checkOverlaps);
  }
  ContainLog->SetVisAttributes(yellow);

    //--Absorber (BeO)
  AbsorberSolid = new G4Tubs("AmBeSolid", 0., AbsorberRadius, 0.5*AbsorberLength, 0., 360.*deg);
  AbsorberLog = new G4LogicalVolume(AbsorberSolid, matAmBe, "AmBeLogical");
  AbsorberLog->SetVisAttributes(red);

    //-- AmBe scorer vertical and perpendicular
  G4VSolid* EnerSphereSolid = new G4Sphere("EnerSphereSolid", WorldSizeXY-1.*mm, WorldSizeXY-0.9*mm, 0., 360.*deg, 0., 180.*deg);
  //G4VSolid* EnerSphereSolid = new G4Sphere("EnerSphereSolid", ContainRadius*3, ContainRadius*3+0.1*mm, 0., 360.*deg, 0., 180.*deg);
  G4LogicalVolume* EnerSphereLogical = new G4LogicalVolume(EnerSphereSolid, matVacuum, "EnerSphereLogical");


  //-- PLACEMENTS --//
  if (fworldType>=1 and fworldType<=3 and fRadioIsotope!="USF" and fRadioIsotope.substr(0,6)!="Single")
  {
    //only source
    new G4PVPlacement(0, G4ThreeVector(), ContainLog, "Container", WorldLog, false, 0, checkOverlaps);
    // Place the absorber inside the inner vacuum when using the X3 casing
    if (fCasingSelection == 1 && X3VacLog)
      new G4PVPlacement(0, G4ThreeVector({0,0,-X3ShiftZ}), AbsorberLog, "AmBe", IT3Log, false, 0, checkOverlaps);
    else
      new G4PVPlacement(0, G4ThreeVector(), AbsorberLog, "AmBe", ContainLog, false, 0, checkOverlaps);
  }

  if (fAzimuthalScoring)//test for neutron spectrum perpendicular and vertical
  {
    if (fworldType>=1 and fworldType<=3 and fRadioIsotope!="USF" and fRadioIsotope.substr(0,6)!="Single")
    {
      new G4PVPlacement(0, G4ThreeVector(), EnerSphereLogical, "EnerSphere", WorldLog, false, 0, checkOverlaps);
    }
  }

  //Show the material table
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  // The Construct() method has to return the final (physical) world volume:
  return worldPhys;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DetectorMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/detector/", "Geometry control");
  fMessenger->DeclareProperty("Debug",fDebug);
  fMessenger->DeclareProperty("RadioIsotope",fRadioIsotope);
  fMessenger->DeclareProperty("AzimuthalScoring",fAzimuthalScoring);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
