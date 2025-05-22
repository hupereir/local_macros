#include <TFile.h>
#include <TH3.h>

#include <tpccalib/TpcSpaceChargeMatrixInversion.h>

#include <algorithm>
#include <utility>
#include <vector>

R__LOAD_LIBRARY(libtpccalib.so)

//_________________________________________
void ExtrapolateDistortionCorrection()
{

  const std::string inputfile = "distortion_maps/average_minus_static_distortion_projected.root";
  const std::string inputfile_cm = "distortion_maps/average_minus_static_distortion_cm.root";
  const std::string outputfile = "distortion_maps/average_minus_static_distortion_extrapolated.root";

  std::cout << "ExtrapolateDistortionCorrection - inputfile: " << inputfile << std::endl;
  std::cout << "ExtrapolateDistortionCorrection - inputfile_cm: " << inputfile_cm << std::endl;
  std::cout << "ExtrapolateDistortionCorrection - outputfile: " << outputfile << std::endl;

  TpcSpaceChargeMatrixInversion matrix_inversion;
  matrix_inversion.load_average_distortion_corrections( inputfile );
  matrix_inversion.load_cm_distortion_corrections( inputfile_cm );
  matrix_inversion.extrapolate_distortion_corrections();
  matrix_inversion.save_distortion_corrections( outputfile );
}
