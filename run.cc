#include "run.hh"

MyRunAction::MyRunAction()
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  //create new output file for every run
  // G4AnalysisManager *man = G4AnalysisManager::Instance();
 
  G4int runID = run->GetRunID();
  count = 0;
  std::stringstream strRunID;
  strRunID << runID;
  
  //man->OpenFile("output"+strRunID.str()+".root");
  //man->CreateNtuple("Counts", "Counts");
  //man->CreateNtupleIColumn("Counts");
  //man->FinishNtuple();
}

void MyRunAction::AddCount()
{
  count += 1;
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
  //write to histogram
  //G4AnalysisManager *man = G4AnalysisManager::Instance();
  //man->Write();
  //man->CloseFile();
  std::cout<<"Count in detector = "<<count<<"\n\n";
}
