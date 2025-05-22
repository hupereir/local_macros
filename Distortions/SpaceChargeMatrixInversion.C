#include <tpccalib/TpcSpaceChargeMatrixInversion.h>

#include <cstdio>
#include <sstream>

R__LOAD_LIBRARY(libtpccalib.so)

//_______________________________________________
// get list of files matching selection
std::vector<std::string> list_files( const std::string& selection )
{
  std::vector<std::string> out;

  std::cout << "list_files - selection: " << selection << std::endl;
  if( selection.empty() ) return out;

  const std::string command = std::string("ls -1 ") + selection;
  auto tmp = popen( command.c_str(), "r" );
  char line[512];
  while( fgets( line, 512, tmp ) )
  {

    std::istringstream istr( line );
    std::string filename;
    istr >> filename;

    if( filename.empty() ) continue;
    if( access( filename.c_str(), R_OK ) ) continue;

    out.push_back( filename );
  }
  pclose( tmp );
  return out;
}

//_______________________________________________
TString SpaceChargeMatrixInversion()
{

  // input files
  // const TString tag = "_flat_acts_truth_notpc_nodistortion";
  // const TString tag = "_flat_acts_truth_notpc_distorted-new";
  const TString tag = "_flat_genfit_truth_notpc_distorted-new";
  const TString inputFile = Form( "DST/CONDOR%s/TpcSpaceChargeMatrices%s_*.root", tag.Data(), tag.Data() );

  // output file
  const TString outputFile = Form( "Rootfiles/Distortions_full%s_mm.root", tag.Data() );
  std::cout << "SpaceChargeMatrixInversion - outputFile: " << outputFile << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion matrix_inversion;

  auto filenames = list_files( inputFile.Data() );
  std::cout << "SpaceChargeMatrixInversion - loaded " << filenames.size() << " files" << std::endl;

  for( const auto& file:filenames )
  { matrix_inversion.add_from_file( file ); }

  matrix_inversion.calculate_distortion_corrections();
  matrix_inversion.save_distortion_corrections( outputFile.Data() );

  return outputFile;

}
