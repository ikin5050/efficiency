#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction, MyRunAction *runAction)
{
  fEventAction = eventAction;
  fRunAction = runAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
  //std::cout<<"in user stepping action\n";
  G4Track *track = step->GetTrack();
  G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());  
  G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
  if(volume == fScoringVolume){
    if(step->IsFirstStepInVolume() && track->GetDynamicParticle()->GetParticleDefinition() == G4Gamma::Definition()){
      //std::cout<<"Inside Scoring Volume"<<"\n\n";
      fRunAction->AddCount();
    }
  }

}
