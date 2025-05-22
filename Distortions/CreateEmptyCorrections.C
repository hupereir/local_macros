#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <g4eval/TrackingEvaluator_hp.h>
#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TTree.h>

#include <Eigen/Dense>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libtpccalib.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

// new bins (from Ross)
constexpr int m_phibins = 80;
constexpr int m_rbins = 52;
constexpr int m_zbins = 160;

// phi range
constexpr float m_phimin = 0;
constexpr float m_phimax = 2.*M_PI;

// r range
constexpr float m_rmin = 20;
constexpr float m_rmax = 78;

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;

//____________________________________
TString CreateEmptyCorrections()
{

  // file definitions
  const auto outputfilename =  "distortion_maps_rec/distortion_corrections_empty.root";
  std::cout << "CreateEmptyCorrections - outputfile: " << outputfilename << std::endl;

  RootFile outputRootFile( outputfilename );

  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = [&]( const TString& name_in, const TString& name_out )
  {
    auto h =  new TH3F( name_in, name_in,
      m_phibins, m_phimin,  m_phimax, 
      m_rbins, m_rmin,  m_rmax, 
      m_zbins, m_zmin,  m_zmax ); 

    h->GetXaxis()->SetTitle( "phi (rad)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "z (cm)" );

    // add source histogram to output
    // outputRootFile.Add( h );

    // split and add
    const auto [hneg, hpos] = TpcSpaceChargeReconstructionHelper::split( h );
    // outputRootFile.Add( hneg );
    // outputRootFile.Add( hpos );

    Utils::PrintAxis( h );

    // copy and add to output
    for( const auto& h: {
      TpcSpaceChargeReconstructionHelper::copy_histogram( h, name_out ),
      TpcSpaceChargeReconstructionHelper::copy_histogram( hneg, Form( "%s_negz", name_out.Data() ) ),
      TpcSpaceChargeReconstructionHelper::copy_histogram( hpos, Form( "%s_posz", name_out.Data() ) ) } )
    {

      Utils::PrintAxis( h );

      // write
      outputRootFile.Add( h );
    }

  };

  std::initializer_list<std::pair<TString, TString>> names =
  {
    // { "hentries_rec", "hentries" },
    { "hDistortionP_rec", "hIntDistortionP" },
    { "hDistortionR_rec", "hIntDistortionR" },
    // { "hDistortionX_rec", "hIntDistortionX" },
    // { "hDistortionY_rec", "hIntDistortionY" },
    { "hDistortionZ_rec", "hIntDistortionZ" }
  };

  for( const auto& [name_in, name_out]:names )
  { process_histogram( name_in, name_out ); }

  return outputfilename;
}
