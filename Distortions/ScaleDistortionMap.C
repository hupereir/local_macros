#include <RootUtil/Utils.h>
#include <RootUtil/RootFile.h>

#include <TH3.h>
#include <TFile.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

TString ScaleDistortionMap()
{
  
  static constexpr double scale = 2.;
   
  const TString inputfile( "distortion_maps/average_minus_static_distortion_converted.root" );
  const TString outputfile( "distortion_maps/average_minus_static_distortion_converted_scaled_x2.root" );
  
  // open input
  auto in = TFile::Open( inputfile );

  // create output rootfile
  RootFile out( outputfile );
  
  for( const auto& hname:{
    "hIntDistortionP_negz","hIntDistortionR_negz","hIntDistortionZ_negz",
    "hIntDistortionP_posz","hIntDistortionR_posz","hIntDistortionZ_posz",
  } )
  {
    auto h = dynamic_cast<TH3*>(in->Get(hname));
    if( !h ) continue;
    std::cout << "ScaleDistortionMap - processing histogram " << hname << std::endl;
    
    h->Scale( scale );

    out.Add( h );
  }
  
  return outputfile;
  
}
