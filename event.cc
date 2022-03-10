#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{
}

MyEventAction::~MyEventAction()
{}

//void MyEventAction::AddCount()
//{
  //count += 1;
//}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
  //std::cout<<"beginning of event action\n";
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  //G4AnalysisManager *man = G4AnalysisManager::Instance();

  //man->FillNtupleIColumn(0, count);
  //man->AddNtupleRow();
}
