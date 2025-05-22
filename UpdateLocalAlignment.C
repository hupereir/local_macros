#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <fstream>
#include <iostream>

R__LOAD_LIBRARY(libtrack_io.so)

//______________________________________________-
class AlignmentParameters
{

  public:

  using list = std::vector<AlignmentParameters>;

  uint32_t hitsetkey = 0;
  double dalpha = 0;
  double dbeta = 0;
  double dgamma = 0;

  double dx = 0;
  double dy = 0;
  double dz = 0;

  double dgalpha = 0;
  double dgbeta = 0;
  double dggamma = 0;

  // add phi rotation (radians)
  void add_global_phi_rotation( double value )
  {

// //  DST/CONDOR_CombinedDataReconstruction_corrected_notpc-rotated2
//     dgalpha += value;
//     dgbeta += value;

    //  DST/CONDOR_CombinedDataReconstruction_corrected_notpc-rotated2
    dgbeta -= value;
    dggamma -= value;

  }

  // streamer
  friend std::istream& operator >> (std::istream& in, AlignmentParameters& param )
  {
    in
      >> param.hitsetkey
      >> param.dalpha >> param.dbeta >> param.dgamma
      >> param.dx >> param.dy >> param.dz
      >> param.dgalpha >> param.dgbeta >> param.dggamma;
    return in;
  }

  // streamer
  friend std::ostream& operator << (std::ostream& out, const AlignmentParameters& param )
  {
    out
      << param.hitsetkey << " "
      << param.dalpha << " " << param.dbeta << " " << param.dgamma << " "
      << param.dx << " " << param.dy << " " << param.dz << " "
      << param.dgalpha << " " << param.dgbeta << " " << param.dggamma;
    return out;
  }

};

//______________________________________________________
void UpdateLocalAlignment()
{

  AlignmentParameters::list parameters;

  const std::string inputfile = "Alignment/localAlignmentParamsFile_orig.txt";
  const std::string outputfile = "Alignment/localAlignmentParamsFile_rotated9.txt";

  std::cout << "UpdateLocalAlignment - inputfile: " << inputfile << std::endl;
  std::cout << "UpdateLocalAlignment - outputfile: " << outputfile << std::endl;

  ifstream in(inputfile.c_str());
  for(;;)
  {
    AlignmentParameters param;
    in >> param;
    if( in )
    {
      parameters.push_back( param );
    } else break;
  }
  std::cout << "ReadLocalAlignment - parameters.size: " << parameters.size() << std::endl;

  // rotate all silicons
  // static const double dphi = -0.03;
  static const double dphi = -0.02;
  for( auto&& parameter:parameters )
  {
    const auto trkrid = TrkrDefs::getTrkrId(parameter.hitsetkey);
    if( trkrid==TrkrDefs::mvtxId || trkrid==TrkrDefs::inttId )
    { parameter.add_global_phi_rotation(dphi); }
  }

  ofstream out(outputfile.c_str());
  for( const auto& parameter:parameters )
  { out << parameter << std::endl; }

  out.close();
}
