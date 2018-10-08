

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "globals.hh"

//------------------------------------------------------------------------------

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim_gen)
:primGenAction(prim_gen)
{ 
  primGenDir = new G4UIdirectory("/primaryGenerator/");
  primGenDir->SetGuidance("UI commands specific to the Primary Generator Action.\n");

  setBeamEnergyCMD = new G4UIcmdWithADoubleAndUnit("/primaryGenerator/setEnergy", this);
  setBeamEnergyCMD->SetGuidance("Set X-ray energy for fixed energy case\n");
  setBeamEnergyCMD->SetGuidance("example: \'/primaryGenerator/setEnergy 3.6 MeV\'\n");
  setBeamEnergyCMD->SetUnitCategory("Energy");
  setBeamEnergyCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  setFieldSizeCMD = new G4UIcmdWithADoubleAndUnit("/primaryGenerator/setFieldSize", this);
  setFieldSizeCMD->SetGuidance("Set X-ray illumination field size\n");
  setFieldSizeCMD->SetGuidance("example: \'/primaryGenerator/setFieldSize 10.44 cm\'\n");
  setFieldSizeCMD->SetUnitCategory("Length");
  setFieldSizeCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  setFixedEnergyFlagCMD = new G4UIcmdWithABool("/primaryGenerator/setFixedEnergyFlag", this);
  setFixedEnergyFlagCMD->SetGuidance("Set fixed energy flag; \'true\' or \'false\'\n");
  setFixedEnergyFlagCMD->SetGuidance("example: \'/primaryGenerator/setFixedEnergyFlag true\' for fixed energy\n");
  setFixedEnergyFlagCMD->SetGuidance("for setting the value of energy see \'/primaryGenerator/setEnergy\'\n");
  setFixedEnergyFlagCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  setParallelBeamFlagCMD = new G4UIcmdWithABool("/primaryGenerator/setParallelBeamFlag", this);
  setParallelBeamFlagCMD->SetGuidance("Set parallel beam flag; \'true\' or \'false\'\n");
  setParallelBeamFlagCMD->SetGuidance("example: \'/primaryGenerator/setParallelBeamFlag true\' for parallel square beam\n");
  setParallelBeamFlagCMD->SetGuidance("for setting the extent/cross-section of the beam see \'/primaryGenerator/setFieldSize\'\n");
  setParallelBeamFlagCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  initializeEnergySpectrumCMD =  new G4UIcmdWithAString("/primaryGenerator/initializeEnergySpectrum", this);
  initializeEnergySpectrumCMD->SetGuidance("Initialize X-ray energy spectrum from a given file path\n");
  initializeEnergySpectrumCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  setBeamPositionOffsetCMD = new G4UIcmdWith3VectorAndUnit("/primaryGenerator/setBeamPositionOffset", this);
  setBeamPositionOffsetCMD->SetGuidance("Set beam origin X Y Z offset;\n");
  setBeamPositionOffsetCMD->SetGuidance("example: \'/primaryGenerator/setBeamPositionOffset x y z unit\'");
  setBeamPositionOffsetCMD->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//------------------------------------------------------------------------------

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete setBeamEnergyCMD;
  delete setFieldSizeCMD;
  delete setFixedEnergyFlagCMD; 
  delete setParallelBeamFlagCMD;
  delete initializeEnergySpectrumCMD;
  delete setBeamPositionOffsetCMD;
  delete primGenDir;
}

//------------------------------------------------------------------------------

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  // get the pointer to the User Interface manager
  //G4UImanager* UI = G4UImanager::GetUIpointer();

  if(command == setBeamEnergyCMD )
    {
      primGenAction->SetParticleEnergy(setBeamEnergyCMD->GetNewDoubleValue(newValue));
      // UI->ApplyCommand("any command here");
    }
  else if ( command == setFieldSizeCMD )
    {
      primGenAction->SetFieldSize(setFieldSizeCMD->GetNewDoubleValue(newValue));
    }
  else if ( command == setFixedEnergyFlagCMD )
    {
      primGenAction->SetFixedEnergyFlag(setFixedEnergyFlagCMD->GetNewBoolValue(newValue));
    }
  else if (command == setParallelBeamFlagCMD )
    {
      primGenAction->SetParallelBeamFlad(setParallelBeamFlagCMD->GetNewBoolValue(newValue));
    }
  else if ( command == initializeEnergySpectrumCMD )
    {
      primGenAction->InitializeSpectrum(newValue);
    }
  else if ( command == setBeamPositionOffsetCMD)
    {
      primGenAction->SetBeamPositionOffset(setBeamPositionOffsetCMD->GetNew3VectorValue(newValue));
    }
}
//------------------------------------------------------------------------------
