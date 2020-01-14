
#include <RootUtil/Stream.h>

#include <fstream>
#include <sstream>

R__LOAD_LIBRARY(libRootUtilBase.so)

void ParseChargedHadronSpectra_2()
{

  std::ifstream in( "data/ChargedHadronSpectra-2.txt" );
  std::string line;

  std::vector<float> pt_vect;
  std::vector<float> yield_vect;
  std::vector<float> stat_vect;

  while( std::getline( in, line ) )
  {
    float pt = 0;
    float yield = 0;
    float stat = 0;

    std::istringstream str( line );
    str >> pt >> yield >> stat;
    if( !str )
    {
      std::cout << "skipping line " << line << std::endl;
      continue;
    }
    pt_vect.push_back( pt );
    yield_vect.push_back( yield );
    stat_vect.push_back( stat );
  }

  // bin limits
  std::vector<float> ptbin_vect = {0.2};
  for( int i = 0; i < pt_vect.size()-1; ++i )
  { ptbin_vect.push_back( 0.5*(pt_vect[i]+pt_vect[i+1] ) ); }

  ptbin_vect.push_back( 3.0 );

  // print
  Stream::PrintVector<float>( "float", "pt", pt_vect, "%.3g" );
  Stream::PrintVector<float>( "float", "ptbins", ptbin_vect, "%.3g" );
  Stream::PrintVector<float>( "float", "yield", yield_vect, "%.3g" );
  Stream::PrintVector<float>( "float", "stat", stat_vect, "%.2g" );
}
