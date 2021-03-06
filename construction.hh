#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume;}

  virtual G4VPhysicalVolume *Construct();

private:
  virtual void ConstructSDandField();

  G4double distance;
  G4double xWorld, yWorld, zWorld, personRadius;

  G4Box *solidWorld;
  G4Box *solidDetector;
  G4LogicalVolume *logicWorld, *logicDetector;
  G4VPhysicalVolume *physWorld, *physDetector;

  G4GenericMessenger *fMessenger;

  G4Material *SiO2, *H2O, *Aerogel, *worldMat, *NaI, *GAGG;
  G4Element *C, *Na, *I, *Gd, *Ce, *O, *Al, *Ga;

  G4LogicalVolume *fScoringVolume;

  void DefineMaterials();
  void ConstructDetector();
};

#endif
