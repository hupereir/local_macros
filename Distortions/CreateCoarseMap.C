#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//___________________________________________________
void FillHistogram( const std::array<TH3*,2>& source, TH3* destination )
{

  // loop over destination bins
  for( int ix = 0; ix < destination->GetNbinsX(); ++ix )
    for( int iy = 0; iy < destination->GetNbinsY(); ++iy )
    for( int iz = 0; iz < destination->GetNbinsZ(); ++iz )
  {
    // get bin center
    const auto x = destination->GetXaxis()->GetBinCenter( ix+1 );
    const auto y = destination->GetYaxis()->GetBinCenter( iy+1 );
    const auto z = std::abs( destination->GetZaxis()->GetBinCenter( iz+1 ) );

    // get content from source histogram
    const auto content = (z<0 ? source[0]:source[1])->Interpolate( x, y, z );
    destination->SetBinContent( ix+1, iy+1, iz+1, content );
  }

}

//___________________________________________________
TString CreateCoarseMap()
{
  
  const TString inputfile( "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps-new/average_minus_static_distortion_converted.root" );
  const TString outputfile( "/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/average_minus_static_distortion_coarse.root" );
  
  std::cout << "CreateCoarseMap - inputfile: " << inputfile << std::endl;
  std::cout << "CreateCoarseMap - outputfile: " << outputfile << std::endl;

  // open input
  auto f = TFile::Open( inputfile );

  // load histograms
  std::array<TH3*, 2> hDPint = {{ dynamic_cast<TH3*>(f->Get("hIntDistortionP_negz")),  dynamic_cast<TH3*>(f->Get("hIntDistortionP_posz")) }};
  std::array<TH3*, 2> hDRint = {{ dynamic_cast<TH3*>(f->Get("hIntDistortionR_negz")),  dynamic_cast<TH3*>(f->Get("hIntDistortionR_posz")) }};
  std::array<TH3*, 2> hDZint = {{ dynamic_cast<TH3*>(f->Get("hIntDistortionZ_negz")),  dynamic_cast<TH3*>(f->Get("hIntDistortionZ_posz")) }};
  Utils::PrintAxis( hDPint[0] );
  Utils::PrintAxis( hDPint[1] );

  // binning and range
  constexpr int m_phibins = 36;
  constexpr int m_rbins = 16;
  constexpr int m_zbins = 80;

  // phi range
  constexpr float m_phimin = 0;
  constexpr float m_phimax = 2.*M_PI;

  // r range
  constexpr float m_rmin = 20;
  constexpr float m_rmax = 78;

  // z range
  constexpr float m_zmin = -105.5;
  constexpr float m_zmax = 105.5;

  // create output histograms
  auto hphi_rec = new TH3F( "hDistortionP_rec", "hDistortionP_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hr_rec = new TH3F( "hDistortionR_rec", "hDistortionR_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hz_rec = new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );

  // copy content
  FillHistogram( hDPint, hphi_rec );
  FillHistogram( hDRint, hr_rec );
  FillHistogram( hDZint, hz_rec );

  RootFile out( outputfile );
  for( const auto& h: {hphi_rec, hr_rec, hz_rec} )
  { out.Add( h ); }

  return outputfile;

}
