#include "generator.hh"


MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particle = particleTable->FindParticle("geantino");

  G4ThreeVector pos(0.,0.,0.);
  G4ThreeVector mom(0.,0.,1.);

  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(0.*GeV);
  fParticleGun->SetParticleDefinition(particle);

  //fGeneralParticleSource = new G4GeneralParticleSource();
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  //delete fGeneralParticleSource;
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  //fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
  G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();

  if(particle==G4Geantino::Geantino())
  {
    G4int Z = 77;
    G4int A = 192;
  
    G4double charge = 0.*eplus;
    G4double energy = 0.*keV;
  
    G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
  
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(charge);
    std::cout<<"\n\n\n\n\nParticle gun setup \n\n\n\n\n\n\n\n\n";
  }
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
