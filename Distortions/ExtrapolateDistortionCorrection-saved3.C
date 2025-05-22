#include <TFile.h>
#include <TH3.h>

#include <algorithm>
#include <utility>
#include <vector>

namespace
{
  /// square
  template<class T>
  inline constexpr T square( T x ) { return x*x; }

  using window_t = std::pair<double, double>;
  using window_list_t = std::vector<window_t>;

  window_t transform_window( const window_t& a )
  { return { a.first+2.*M_PI, a.second+2.*M_PI }; }

  // detector geometry
  const window_t phi_window_central=transform_window({-1.742,-1.43979});
  const window_t phi_window_east=transform_window({-2.27002,-1.9673});
  const window_t phi_window_west=transform_window({-1.21452,-0.911172});

  const window_list_t theta_window_central={{-0.918257,-0.613136}, {-0.567549,-0.031022}, {0.0332154,0.570419}, {0.613631,0.919122}};
  const window_list_t theta_window_east={{-0.636926,-0.133603}, {0.140678,0.642714}};
  const window_list_t theta_window_west={{-0.643676,-0.141004}, {0.13485,0.640695}};

  class range_ftor_t
  {
    public:
    range_ftor_t( double value ): m_value( value ){};
    bool operator ()( const window_t& range )
    { return m_value > range.first && m_value < range.second; }

    private:
    double m_value = 0;
  };

  bool in_range( const double& value, const window_t range )
  { return range_ftor_t(value)(range); }

  bool in_range( const double& value, const window_list_t range_list )
  { return std::any_of( range_list.begin(), range_list.end(), range_ftor_t(value)); }
}

//_________________________________________
/* create a mask histogram that corresponds to TPOT acceptance, using the bins from the reference histogram */
TH3* create_mask( TH3* source )
{
  auto hmask = static_cast<TH3*>( source->Clone( "mask" ) );
  hmask->Reset();

  // loop over bins
  for( int ip = 0; ip < source->GetNbinsX(); ++ip )
    for( int ir = 0; ir < source->GetNbinsY(); ++ir )
    for( int iz = 0; iz < source->GetNbinsZ(); ++iz )
  {
    const double phi = source->GetXaxis()->GetBinCenter(ip+1);
    const double r = source->GetYaxis()->GetBinCenter(ir+1);
    const double z = source->GetZaxis()->GetBinCenter(iz+1);
    const double theta = std::atan2( z, r );
    if( in_range( phi, phi_window_central ) && in_range( theta, theta_window_central )  )
    { hmask->SetBinContent( ip+1, ir+1, iz+1, 1 ); }
  }

  return hmask;
}

//_________________________________________
void extrapolate_phi( TH3* source, TH3* destination, TH2* h_cm, TH3* mask )
{

  // loop over phi bins
  for( int ip = 0; ip < source->GetNbinsX(); ++ip )
    for( int ir = 0; ir < source->GetNbinsY(); ++ir )
    for( int iz = 0; iz < source->GetNbinsZ(); ++iz )
  {
    bool in_range = mask->GetBinContent( ip+1, ir+1, iz+1 ) > 0;
    if( in_range )
    {
      destination->SetBinContent( ip+1, ir+1, iz+1, source->GetBinContent(ip+1, ir+1, iz+1) );
      continue;
    }

    // phi not in TPOT range, rotate by steps of 2pi/12 (= one TPC sector) until found in range
    const double phi = source->GetXaxis()->GetBinCenter(ip+1);
    static constexpr int n_sectors = 12;
    for( int sector = 1; sector < n_sectors; ++sector )
    {
      // get ref phi
      double phi_ref = phi + 2.*M_PI*sector/n_sectors;
      while( phi_ref >= 2*M_PI ) { phi_ref -= 2*M_PI; };

      int ip_ref = source->GetXaxis()->FindBin( phi_ref ) - 1;
      in_range = mask->GetBinContent( ip_ref+1, ir+1, iz+1 ) > 0;
      if( !in_range ) continue;

      // get normalization factor from CM histograms
      double scale = 1;
      if( h_cm )
      {
        const double distortion_local = h_cm->GetBinContent( ip+1, ir+1 );
        const double distortion_ref = h_cm->GetBinContent( ip_ref+1, ir+1 );
        scale = distortion_local/distortion_ref;
      }

      destination->SetBinContent( ip+1, ir+1, iz+1, scale*source->GetBinContent( ip_ref+1, ir+1, iz+1 ) );
    }
  }
}

//_________________________________________
void extrapolate_phi2( TH3* source, TH3* destination, TH3* mask )
{

  // loop over phi bins
  for( int ir = 0; ir < source->GetNbinsY(); ++ir )
    for( int iz = 0; iz < source->GetNbinsZ(); ++iz )
  {
    int ip_min = -1;
    for( int ip = 0; ip < source->GetNbinsX(); ++ip )
    {
      bool in_range = mask->GetBinContent( ip+1, ir+1, iz+1 ) > 0;
      if( in_range )
      {
        destination->SetBinContent( ip+1, ir+1, iz+1, source->GetBinContent(ip+1, ir+1, iz+1) );
        ip_min = ip;
        continue;
      }

      // ip is not in range. Check if ip_min was set. If not, skip this bin. If yes, find next bin in range
      if( ip_min < 0 ) continue;

      int ip_max = ip+1;
      for( ; ip_max < source->GetNbinsX(); ++ip_max )
      {
        in_range = mask->GetBinContent( ip_max+1, ir+1, iz+1 ) > 0;
        if( in_range ) break;
      }

      // check that a valid bin was found
      if( ip_max == source->GetNbinsX() ) continue;

      // do the interpolation
      const double phi_min =  source->GetXaxis()->GetBinCenter(ip_min+1);
      const double phi_max =  source->GetXaxis()->GetBinCenter(ip_max+1);
      const double phi = source->GetXaxis()->GetBinCenter(ip+1);

      const double distortion_min = source->GetBinContent( ip_min+1, ir+1, iz+1 );
      const double distortion_max = source->GetBinContent( ip_max+1, ir+1, iz+1 );
      const double distortion = distortion_min + (phi-phi_min)*(distortion_max - distortion_min)/(phi_max-phi_min);

      destination->SetBinContent( ip+1, ir+1, iz+1, distortion );
    }
  }
}

//_________________________________________
void extrapolate_z( TH3* source, TH3* destination, TH3* mask )
{
  // loop over phi bins
  for( int ir = 0; ir < source->GetNbinsY(); ++ir )
    for( int ip = 0; ip < source->GetNbinsX(); ++ip )
  {
    int iz_min = -1;
    for( int iz = 0; iz < source->GetNbinsZ(); ++iz )
    {
      bool in_range = mask->GetBinContent( ip+1, ir+1, iz+1 ) > 0;
      if( in_range )
      {
        destination->SetBinContent( ip+1, ir+1, iz+1, source->GetBinContent(ip+1, ir+1, iz+1) );
        iz_min = iz;
        continue;
      }

      // iz is not in range. Check if iz_min was set. If not, skip this bin. If yes, find next bin in range
      if( iz_min < 0 ) continue;

      int iz_max = iz+1;
      for( ; iz_max < source->GetNbinsZ(); ++iz_max )
      {
        in_range = mask->GetBinContent( ip+1, ir+1, iz_max+1 ) > 0;
        if( in_range ) break;
      }

      // check interpolation range validity
      if( iz_max == source->GetNbinsZ() ) continue;

      // do the interpolation
      const double z_min =  source->GetZaxis()->GetBinCenter(iz_min+1);
      const double z_max =  source->GetZaxis()->GetBinCenter(iz_max+1);
      const double z = source->GetZaxis()->GetBinCenter(iz+1);
      const double alpha_min = (z_max-z)/(z_max-z_min);
      const double alpha_max = (z-z_min)/(z_max-z_min);

      const double distortion_min = source->GetBinContent( ip+1, ir+1, iz_min+1 );
      const double distortion_max = source->GetBinContent( ip+1, ir+1, iz_max+1 );
      const double distortion = alpha_min*distortion_min + alpha_max*distortion_max;

      const double error_min = source->GetBinError( ip+1, ir+1, iz_min+1 );
      const double error_max = source->GetBinError( ip+1, ir+1, iz_max+1 );
      const double error = std::sqrt(square(alpha_min*error_min) + square(alpha_max*error_max));

      destination->SetBinContent( ip+1, ir+1, iz+1, distortion );
      destination->SetBinError( ip+1, ir+1, iz+1, error );
    }
  }
}

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
  auto hentries_posz = create_mask( static_cast<TH3*>( in->Get("hIntDistortionR_posz" ) ) );
  hentries_posz->SetName( "hentries_posz" );

  auto hentries_posz_extrapz = static_cast<TH3*>( hentries_posz->Clone() );
  extrapolate_z( hentries_posz, hentries_posz_extrapz, hentries_posz );

  auto hentries_posz_extrap = static_cast<TH3*>( hentries_posz->Clone() );
  extrapolate_phi( hentries_posz_extrapz, hentries_posz_extrap, nullptr, hentries_posz_extrapz );

  for( const TString& hname: {"hIntDistortionR_posz", "hIntDistortionP_posz", "hIntDistortionZ_posz" } )
  {
    auto h = static_cast<TH3*>( in->Get(hname) );
    auto h_cm = static_cast<TH2*>( in_cm->Get(hname) );
    auto h_extrap = static_cast<TH3*>( h->Clone( TString( h->GetName() )+"_copy" ) );
    extrapolate_z( h, h_extrap, hentries_posz );
    extrapolate_phi( h_extrap, h_extrap, h_cm, hentries_posz_extrapz );
    extrapolate_phi2( h_extrap, h_extrap, hentries_posz_extrap );

    out->cd();
    h_extrap->Write( hname );
  }

  // negative z
  auto hentries_negz = create_mask( static_cast<TH3*>( in->Get("hIntDistortionR_negz" ) ) );
  hentries_negz->SetName( "hentries_negz" );

  auto hentries_negz_extrapz = static_cast<TH3*>(hentries_negz->Clone());
  extrapolate_z( hentries_negz, hentries_negz_extrapz, hentries_negz );

  auto hentries_negz_extrap = static_cast<TH3*>(hentries_negz->Clone());
  extrapolate_phi( hentries_negz_extrapz, hentries_negz_extrap, nullptr, hentries_negz_extrapz );

  for( const TString& hname: {"hIntDistortionR_negz", "hIntDistortionP_negz", "hIntDistortionZ_negz" } )
  {
    auto h = static_cast<TH3*>( in->Get(hname) );
    auto h_cm = static_cast<TH2*>( in_cm->Get(hname) );
    auto h_extrap = static_cast<TH3*>( h->Clone( TString( h->GetName() )+"_copy" ) );
    extrapolate_z( h, h_extrap, hentries_negz );
    extrapolate_phi( h_extrap, h_extrap, h_cm, hentries_negz_extrapz );
    extrapolate_phi2( h_extrap, h_extrap, hentries_negz_extrap );

    out->cd();
    h_extrap->Write( hname );
  }

  // also write extrapolated masks
  out->cd();
//   hentries_posz_extrap->Write("hentries_posz");
//   hentries_negz_extrap->Write("hentries_negz");
  out->Close();

}
