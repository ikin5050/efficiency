#include "construction.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iterator>

template <typename A, typename B>
void zip(
    const std::vector<A> &a, 
    const std::vector<B> &b, 
    std::vector<std::pair<A,B>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

// Write the first and second element of the pairs in 
// the given zipped vector into a and b. (This assumes 
// that the vectors have equal length)
template <typename A, typename B>
void unzip(
    const std::vector<std::pair<A, B>> &zipped, 
    std::vector<A> &a, 
    std::vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}



template <class Container>
void Split(const std::string& str, Container& cont)
{
  std::istringstream iss(str);
  std::copy(std::istream_iterator<std::string>(iss),
	    std::istream_iterator<std::string>(),
	    std::back_inserter(cont));
}


void GetOP(std::string fname, std::vector<G4double> *WL, std::vector<G4double> *Em, std::vector<double> *Ridx){
  double c_spd = 299792458;
  double h_planck =  4.135667516E-6; //ev for nm
  std::fstream fin(fname.c_str());
  std::string line;
  while(std::getline(fin, line)){
      std::vector<std::string> words;
      Split(line,words);
      WL->push_back(c_spd*h_Planck/std::stod(words.at(0)));
      Em->push_back(std::stod(words.at(1)));
      Ridx->push_back(1.9);
    }
  

  fin.close();
  

}

MyDetectorConstruction::MyDetectorConstruction()
{
  fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
  
  fMessenger->DeclareProperty("distance", distance, "distance from origin");
  
  distance = 1.;
  
  DefineMaterials();
  
  xWorld = 0.5*m;
  yWorld = 0.5*m;
  zWorld = 0.5*m;
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

void MyDetectorConstruction::DefineMaterials()
{ 
  G4NistManager *nist = G4NistManager::Instance();

  SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
  SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
  SiO2->AddElement(nist->FindOrBuildElement("O"), 2);

  H2O = new G4Material("H20", 1*g/cm3, 2);
  H2O->AddElement(nist->FindOrBuildElement("H"), 2);
  H2O->AddElement(nist->FindOrBuildElement("O"), 1);

  C = nist->FindOrBuildElement("C");

  Na = nist->FindOrBuildElement("Na");
  I = nist->FindOrBuildElement("I");
  NaI = new G4Material("NaI", 3.67*g/cm3, 2);
  NaI->AddElement(Na, 1);
  NaI->AddElement(I, 1);

 
  Gd = nist->FindOrBuildElement("Gd");
  Al = nist->FindOrBuildElement("Al");
  Ga = nist->FindOrBuildElement("Ga");
  O = nist->FindOrBuildElement("O");
  Ce = nist->FindOrBuildElement("Ce");
  GAGG = new G4Material("GAGG", 6.63*g/cm3, 5);
  GAGG->AddElement(Gd, 50.3918*perCent);
  GAGG->AddElement(Al, 5.761*perCent);
  GAGG->AddElement(Ga, 22.3443*perCent);
  GAGG->AddElement(O, 20.5029*perCent);
  GAGG->AddElement(Ce, 0.01);	   

  std::vector<G4double> GAGG_Wavelength;
  std::vector<G4double> GAGG_Emission;
  std::vector<G4double> GAGG_Rindex;
  
  std::string fname = "./OpticalProperties/GAGG_Emission.csv";
  
  GetOP(fname,&GAGG_Wavelength,&GAGG_Emission,&GAGG_Rindex);
  
  std::vector<std::pair<G4double,G4double>> zipped;
  
  zip(GAGG_Wavelength,GAGG_Emission,zipped);
  
  std::sort(std::begin(zipped), std::end(zipped), 
	     [&](const auto& a, const auto& b)
	     {
  	      return a.first > b.first;
  	    });
   
   unzip(zipped, GAGG_Wavelength,GAGG_Emission); //put things in order
   
   G4MaterialPropertiesTable* MPT_GAGG = new G4MaterialPropertiesTable();
   MPT_GAGG->AddProperty("SLOWCOMPONENT", GAGG_Wavelength,GAGG_Emission)->SetSpline(true);
   MPT_GAGG->AddProperty("RINDEX", GAGG_Wavelength,GAGG_Rindex)->SetSpline(true);
   MPT_GAGG->AddConstProperty("SCINTILLATIONYIELD", 50000. / MeV);
   MPT_GAGG->AddConstProperty("RESOLUTIONSCALE", 1.0);
   MPT_GAGG->AddConstProperty("SLOWTIMECONSTANT", 88. * ns);
   
   GAGG->SetMaterialPropertiesTable(MPT_GAGG);
   //std::cout<<"\n\n\n\n\nGAGG Material Table Set\n\n\n\n\n";

   worldMat = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
}


void MyDetectorConstruction::ConstructDetector()
{
  solidDetector = new G4Box("solidRadiator", 1*cm, 1*cm, 1*cm);
  logicDetector = new G4LogicalVolume(solidDetector, GAGG, "logicDetector");
  fScoringVolume = logicDetector;
  physDetector = new G4PVPlacement(0, G4ThreeVector(distance*cm,0.,0.), logicDetector, "physDetector", logicWorld, false, 0, true);
}


G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);

  logicWorld = new G4LogicalVolume(solidWorld, worldMat,"logicWorld");
  
  physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

  ConstructDetector();
  std::cout<<"\n\nDistance of detector from origin: "<<distance<<"cm \n ";
  
  return physWorld;
  
}


void MyDetectorConstruction::ConstructSDandField()
{
  MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");
  logicDetector->SetSensitiveDetector(sensDet);
  std::cout<<"\n\n\n\nSensitive detector constructed  \n\n\n\n";
}
