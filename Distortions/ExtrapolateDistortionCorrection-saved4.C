#include <TFile.h>
#include <TH3.h>

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <algorithm>
#include <utility>
#include <vector>

R__LOAD_LIBRARY(libtpccalib.so)

//_________________________________________
void ExtrapolateDistortionCorrection()
{

  const TString inputfile = "distortion_maps/average_minus_static_distortion_projected.root";
  const TString inputfile_cm = "distortion_maps/average_minus_static_distortion_cm.root";
  const TString outputfile = "distortion_maps/average_minus_static_distortion_extrapolated.root";

  std::cout << "ExtrapolateDistortionCorrection - inputfile: " << inputfile << std::endl;
  std::cout << "ExtrapolateDistortionCorrection - inputfile_cm: " << inputfile_cm << std::endl;
  std::cout << "ExtrapolateDistortionCorrection - outputfile: " << outputfile << std::endl;

  // open TFiles
  auto in = TFile::Open( inputfile );
  auto in_cm = TFile::Open( inputfile_cm );
  auto out = TFile::Open( outputfile, "RECREATE" );

  // positive z
  auto hmask_posz = static_cast<TH3*>( in->Get("hIntDistortionR_posz" )->Clone( "hentries_posz" ) );
  TpcSpaceChargeReconstructionHelper::create_tpot_mask( hmask_posz );

  auto hmask_posz_extrap_z = static_cast<TH3*>( hmask_posz->Clone( "hentries_posz_extrap_z" ) );
  TpcSpaceChargeReconstructionHelper::extrapolate_z( hmask_posz_extrap_z, hmask_posz );

  auto hmask_posz_extrap_p = static_cast<TH3*>( hmask_posz_extrap_z->Clone( "hentries_posz_extrap_p" ) );
  TpcSpaceChargeReconstructionHelper::extrapolate_phi1( hmask_posz_extrap_p, nullptr, hmask_posz_extrap_z );

  for( const TString& hname: {"hIntDistortionR_posz", "hIntDistortionP_posz", "hIntDistortionZ_posz" } )
  {
    auto h = static_cast<TH3*>( in->Get(hname) );
    auto h_cm = static_cast<TH2*>( in_cm->Get(hname) );
    auto h_extrap = static_cast<TH3*>( h->Clone( TString( h->GetName() )+"_extrap" ) );
    TpcSpaceChargeReconstructionHelper::extrapolate_z( h_extrap, hmask_posz );
    TpcSpaceChargeReconstructionHelper::extrapolate_phi1( h_extrap, h_cm, hmask_posz_extrap_z );
    TpcSpaceChargeReconstructionHelper::extrapolate_phi2( h_extrap, hmask_posz_extrap_p );

    out->cd();
    h_extrap->Write( hname );
  }

  // negative z
  auto hmask_negz = static_cast<TH3*>( in->Get("hIntDistortionR_negz" )->Clone( "hentries_negz" ) );
  TpcSpaceChargeReconstructionHelper::create_tpot_mask( hmask_negz );

  auto hmask_negz_extrap_z = static_cast<TH3*>( hmask_negz->Clone( "hentries_negz_extrap_z" ) );
  TpcSpaceChargeReconstructionHelper::extrapolate_z( hmask_negz_extrap_z, hmask_negz );

  auto hmask_negz_extrap_p = static_cast<TH3*>( hmask_negz_extrap_z->Clone( "hentries_negz_extrap_p" ) );
  TpcSpaceChargeReconstructionHelper::extrapolate_phi1( hmask_negz_extrap_p, nullptr, hmask_negz_extrap_z );

  for( const TString& hname: {"hIntDistortionR_negz", "hIntDistortionP_negz", "hIntDistortionZ_negz" } )
  {
    auto h = static_cast<TH3*>( in->Get(hname) );
    auto h_cm = static_cast<TH2*>( in_cm->Get(hname) );
    auto h_extrap = static_cast<TH3*>( h->Clone( TString( h->GetName() )+"_extrap" ) );
    TpcSpaceChargeReconstructionHelper::extrapolate_z( h_extrap, hmask_negz );
    TpcSpaceChargeReconstructionHelper::extrapolate_phi1( h_extrap, h_cm, hmask_negz_extrap_z );
    TpcSpaceChargeReconstructionHelper::extrapolate_phi2( h_extrap, hmask_negz_extrap_p );

    out->cd();
    h_extrap->Write( hname );
  }

  // also write extrapolated masks
  out->cd();
//   hmask_posz_extrap->Write("hentries_posz");
//   hmask_negz_extrap->Write("hentries_negz");
  out->Close();

}
