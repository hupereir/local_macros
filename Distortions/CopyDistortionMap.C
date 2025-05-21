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

//____________________________________
bool symetrize_histogram( TH3* h, int phibin )
{
  const auto phibins = h->GetXaxis()->GetNbins();
  const auto rbins = h->GetYaxis()->GetNbins();
  const auto zbins = h->GetZaxis()->GetNbins();

  // fill center
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {
    if( iphi+1 != phibin )
    {
      h->SetBinContent( iphi+1, ir+1, iz+1, h->GetBinContent( phibin, ir+1, iz+1 ) );
      h->SetBinError( iphi+1, ir+1, iz+1, h->GetBinError( phibin, ir+1, iz+1 ) );
    }
  }

  return true;
}

//____________________________________
TString CopyDistortionMap()
{

  
  // TString tag = "_acts_truth_notpc_nodistortion";
  // TString tag = "_acts_truth_notpc_distorted";
  TString tag = "_genfit_truth_notpc_nodistortion";
  // TString tag = "_genfit_truth_notpc_distorted";
  const auto inputfilename = Form( "Rootfiles/Distortions_full_realistic_micromegas_mm%s.root", tag.Data() );
  const auto outputfilename = Form( "distortion_maps-new/Distortions_full_realistic_micromegas_mm%s.root", tag.Data() ); 

  std::cout << "CopyDistortionMap - inputfile: " << inputfilename << std::endl;
  std::cout << "CopyDistortionMap - outputfile: " << outputfilename << std::endl;

  RootFile outputRootFile( outputfilename );

  auto in = TFile::Open( inputfilename );

  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = [&]( const TString& name_in, const TString& name_out )
  {
    auto h =  static_cast<TH3*>( in->Get( name_in ) );
    if( !h )
    {
      std::cout << "CopyDistortionMap - cannot find histogram named: " << name_in << std::endl;
      return;
    } else {
      std::cout << "CopyDistortionMap - processing: " << name_in << std::endl;
    }

    // add source histogram to output
    // outputRootFile.Add( h );

    // split and add
    const auto [hneg, hpos] = TpcSpaceChargeReconstructionHelper::split( h );
    // outputRootFile.Add( hneg );
    // outputRootFile.Add( hpos );

    Utils::PrintAxis( h );

    // copy and add to output
    for( const auto& h: {
      // TpcSpaceChargeReconstructionHelper::copy_histogram( h, name_out ),
      TpcSpaceChargeReconstructionHelper::copy_histogram( hneg, Form( "%s_negz", name_out.Data() ) ),
      TpcSpaceChargeReconstructionHelper::copy_histogram( hpos, Form( "%s_posz", name_out.Data() ) ) } )
    {

      Utils::PrintAxis( h );

//       // symetrize along phi
//       if( false )
//       {
//         static constexpr int phibin_rec = 11;
//         symetrize_histogram( h, phibin_rec+1 );
//       }

      // write
      outputRootFile.Add( h );
    }

  };

  std::initializer_list<std::pair<TString, TString>> names =
  {
    { "hentries_rec", "hentries" },
    { "hDistortionP_rec", "hIntDistortionP" },
    { "hDistortionR_rec", "hIntDistortionR" },
//     { "hDistortionX_rec", "hIntDistortionX" },
//     { "hDistortionY_rec", "hIntDistortionY" },
    { "hDistortionZ_rec", "hIntDistortionZ" }
  };

  for( const auto& [name_in, name_out]:names )
  { process_histogram( name_in, name_out ); }

  return outputfilename;
}
