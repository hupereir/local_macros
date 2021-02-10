#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

#include <Eigen/Dense>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

//____________________________________
TH3* create_histogram( TH3* hin, const TString& name )
{
  std::array<int, 3> bins;
  std::array<double, 3> x_min;
  std::array<double, 3> x_max;

  int index = 0;
  for( const auto axis:{ hin->GetXaxis(), hin->GetYaxis(), hin->GetZaxis() } )
  {
    const auto bin_width = (axis->GetXmax() - axis->GetXmin())/axis->GetNbins();
    bins[index] = axis->GetNbins()+2;
    x_min[index] = axis->GetXmin()-bin_width;
    x_max[index] = axis->GetXmax()+bin_width;
    ++index;
  }

  auto hout = new TH3F( name, name,
    bins[0], x_min[0], x_max[0],
    bins[1], x_min[1], x_max[1],
    bins[2], x_min[2], x_max[2] );

  hout->GetXaxis()->SetTitle( "#phi (rad)" );
  hout->GetYaxis()->SetTitle( "r (cm)" );
  hout->GetZaxis()->SetTitle( "z (cm)" );

  return hout;

}

//____________________________________
bool copy_histogram_content( TH3* hin, TH3* hout )
{
  if( !(hin && hout ) ) return false;

  const auto phibins = hin->GetXaxis()->GetNbins();
  const auto rbins = hin->GetYaxis()->GetNbins();
  const auto zbins = hin->GetZaxis()->GetNbins();

  // fill center
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {
    hout->SetBinContent( iphi+2, ir+2, iz+2, hin->GetBinContent( iphi+1, ir+1, iz+1 ) );
    hout->SetBinError( iphi+2, ir+2, iz+2, hin->GetBinError( iphi+1, ir+1, iz+1 ) );
  }

  // fill guarding phi bins
  for( int ir = 0; ir < rbins+2; ++ir )
    for( int iz = 0; iz < zbins+2; ++iz )
  {
    hout->SetBinContent( 1, ir+1, iz+1, hout->GetBinContent( 2, ir+1, iz+1 ) );
    hout->SetBinError( 1, ir+1, iz+1, hout->GetBinError( 2, ir+1, iz+1 ) );

    hout->SetBinContent( phibins+2, ir+1, iz+1, hout->GetBinContent( phibins+1, ir+1, iz+1 ) );
    hout->SetBinError( phibins+2, ir+1, iz+1, hout->GetBinError( phibins+1, ir+1, iz+1 ) );
  }

  // fill guarding r bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int iz = 0; iz < zbins+2; ++iz )
  {
    hout->SetBinContent( iphi+1, 1, iz+1, hout->GetBinContent( iphi+1, 2, iz+1 ) );
    hout->SetBinError( iphi+1, 1, iz+1, hout->GetBinError( iphi+1, 2, iz+1 ) );

    hout->SetBinContent( iphi+1, rbins+2, iz+1, hout->GetBinContent( iphi+1, rbins+1, iz+1 ) );
    hout->SetBinError( iphi+1, rbins+1, iz+1, hout->GetBinError( iphi+1, rbins+1, iz+1 ) );
  }

  // fill guarding z bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int ir = 0; ir < rbins+2; ++ir )
  {
    hout->SetBinContent( iphi+1, ir+1, 1, hout->GetBinContent( iphi+1, ir+1, 2 ) );
    hout->SetBinError( iphi+1, ir+1, 1, hout->GetBinError( iphi+1, ir+1, 2 ) );

    hout->SetBinContent( iphi+1, ir+1, zbins+2, hout->GetBinContent( iphi+1, ir+1, zbins+1 ) );
    hout->SetBinError( iphi+1, ir+1, zbins+2, hout->GetBinError( iphi+1, ir+1, zbins+1 ) );
  }

  return true;

}

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

  // file definitions
  // const TString tag = "_full_realistic_micromegas_mm-new-loose_extrapolated";
  // const TString tag = "_full_realistic_micromegas_mm-new_extrapolated";
  const TString tag = "_full_realistic_micromegas_mm-coarse-new2_extrapolated";

  const auto inputfilename = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  const auto outputfilename =  Form( "distortion_maps_rec/Distortions%s.root", tag.Data() );

//   const auto inputfilename = "Rootfiles/fluct_average-coarse.root";
//   const auto outputfilename = "distortion_maps/fluct_average-coarse.root";

//   const auto inputfilename = "Rootfiles/empty.root";
//   const auto outputfilename = "distortion_maps/empty.root";

  std::cout << "CopyDistortionMap - inputfile: " << inputfilename << std::endl;
  std::cout << "CopyDistortionMap - outputfile: " << outputfilename << std::endl;

  RootFile outputRootFile( outputfilename );

  auto in = TFile::Open( inputfilename );

  std::initializer_list<std::pair<TString, TString>> names = {
    { "hentries_rec", "hentries" },
    { "hDistortionP_rec", "hIntDistortionP" },
    { "hDistortionR_rec", "hIntDistortionR" },
    { "hDistortionX_rec", "hIntDistortionX" },
    { "hDistortionY_rec", "hIntDistortionY" },
    { "hDistortionZ_rec", "hIntDistortionZ" }
  };

  // copy source histogram
  for( const auto& pair:names )
  {
    auto hin =  static_cast<TH3*>( in->Get( pair.first ) );
    outputRootFile.Add( hin );
  }
    
  // generate new histograms and copy
  for( const auto& pair:names )
  {
    auto hin =  static_cast<TH3*>( in->Get( pair.first ) );
    if( !hin )
    {
      std::cout << "CopyDistortionMap - cannot find histogram named: " << pair.first << std::endl;
      continue;
    } else std::cout << "CopyDistortionMap - processing: " << pair.first << std::endl;

    auto hout = create_histogram( hin, pair.second );
    copy_histogram_content( hin, hout );

    // symetrize along phi
    if( false )
    {
      static constexpr int phibin_rec = 11;
      symetrize_histogram( hout, phibin_rec+1 );
    }

    outputRootFile.Add( hout );

    Utils::PrintAxis( hin );
    Utils::PrintAxis( hout );

  }

  // also copy entries
  auto hentries = dynamic_cast<TH3*>(in->Get("hentries_rec"));
  outputRootFile.Add(hentries );

  return outputfilename;
}
