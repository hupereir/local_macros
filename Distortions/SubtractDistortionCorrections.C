#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

TString SubtractDistortionCorrections()
{
  
  const TString firstInputFile( "Rootfiles/Distortions_full_realistic_micromegas_mm-coarse.root" );
  const TString secondInputFile( "Rootfiles/Distortions_full_realistic_micromegas_mm-empty.root" );
  const TString outputFile( "Rootfiles/Distortions_full_realistic_micromegas_mm-coarse-subtracted.root" );
  
  std::cout << "SubtractDistortionCorrections - firstInputFile: " <<  firstInputFile << std::endl;
  std::cout << "SubtractDistortionCorrections - secondInputFile: " <<  secondInputFile << std::endl;
  std::cout << "SubtractDistortionCorrections - outputFile: " <<  outputFile << std::endl;
  
  auto firstTFile = TFile::Open( firstInputFile );
  auto secondTFile = TFile::Open( secondInputFile );
  
  RootFile outputTFile( outputFile );
  
  // read entries from first file and add to output
  outputTFile.Add( static_cast<TH3*>( firstTFile->Get( "hentries_rec" ) ) );
  
  
  for( const auto name: { "hDistortionP_rec", "hDistortionR_rec", "hDistortionZ_rec" } )
  {
    auto hfirst = static_cast<TH3*>( firstTFile->Get( name ) );
    auto hsecond = static_cast<TH3*>( secondTFile->Get( name ) );
    hfirst->Add( hsecond, -1 );
    outputTFile.Add( hfirst );    
  }  
  
  return outputFile;
  
}
