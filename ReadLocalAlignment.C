#include <TrackBase/TrkrDefs.h>

#include <cmath>
#include <fstream>
#include <iostream>

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

  // add phi rotation (radians)
  bool add_phi_rotation( double value )
  {
    dalpha += value;
    dbeta += value-M_PI/2;
  }

  // streamer
  friend std::istream& operator >> (std::istream& in, AlignmentParameters& param )
  {
    double dummy;
    in
      >> param.hitsetkey
      >> param.dalpha >> param.dbeta >> param.dgamma
      >> param.dx >> param.dy >> param.dz
      >> dummy >> dummy >> dummy;
    return in;
  }

  // streamer
  friend std::ostream& operator << (std::ostream& out, const AlignmentParameters& param )
  {
    const double dummy = 0;
    out
      << param.hitsetkey << " "
      << param.dalpha << " " << param.dbeta << " " << param.dgamma << " "
      << param.dx << " " << param.dy << " " << param.dz << " "
      << dummy << " " << dummy << " " << dummy;
    return out;
  }

};

//______________________________________________________
void ReadLocalAlignment()
{
  const std::string inputfile = "localAlignmentParamsFile.txt";
  const std::string outputfile = "localAlignmentParamsFile_new.txt";

  // read parameters from input file
  AlignmentParameters::list parameters;

  {
    ifstream in( inputfile.c_str() );
    for(;;)
    {
      AlignmentParameters param;
      in >> param;
      if( in )
      {
        parameters.push_back( param );
      } else break;
    }
  }

  std::cout << "ReadLocalAlignment - parameters.size: " << parameters.size() << std::endl;
}
