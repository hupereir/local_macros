#include <TFile.h>
#include <TH3.h>

#include <algorithm>
#include <utility>
#include <vector>

namespace
{
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
void ExtrapolateDistortionCorrection_saved()
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

  // loop over histograms
  for( const TString& hname: {
    "hentries_posz", "hentries_negz",
    "hIntDistortionR_posz",  "hIntDistortionR_negz",
    "hIntDistortionP_posz",  "hIntDistortionP_negz",
    "hIntDistortionZ_posz",  "hIntDistortionZ_negz" } )
  {
    auto h = static_cast<TH3*>( in->Get(hname) );
    auto h_cm = static_cast<TH2*>( in_cm->Get(hname) );

    // clone
    auto h_copy = static_cast<TH3*>(h->Clone(hname+"_copy"));
    h_copy->Reset();

    // loop over phi bins
    for( int ip = 0; ip < h->GetNbinsX(); ++ip )
    {
      const double phi = h->GetXaxis()->GetBinCenter(ip+1);

      // phi is in TPOT range. Copy all r,z bins unchanged
      if( in_range( phi, phi_window_central ) )
      {
        for( int ir = 0; ir < h->GetNbinsY(); ++ir )
          for( int iz = 0; iz < h->GetNbinsZ(); ++iz )
        { h_copy->SetBinContent( ip+1, ir+1, iz+1, h->GetBinContent(ip+1, ir+1, iz+1) ); }
        continue;
      }

      // phi not in TPOT range, rotate by steps of 2pi/12 (= one TPC sector) until found in range
      static constexpr int n_sectors = 12;
      for( int sector = 1; sector < n_sectors; ++sector )
      {
        // get ref phi
        double phi_ref = phi + 2.*M_PI*sector/n_sectors;
        while( phi_ref >= 2*M_PI ) { phi_ref -= 2*M_PI; };

        if( !in_range( phi_ref, phi_window_central ) ) continue;

        // get bin number (counting from zero) mathing found sector
        int ip_ref = h->GetXaxis()->FindBin( phi_ref ) - 1;

        // find destination sector, loop over r and z bins
        for( int ir = 0; ir < h->GetNbinsY(); ++ir )
        {
          // get normalization factor from CM histograms
          double scale = 1;
          if( h_cm )
          {
            const double distortion_local = h_cm->GetBinContent( ip+1, ir+1 );
            const double distortion_ref = h_cm->GetBinContent( ip_ref+1, ir+1 );
            scale = distortion_local/distortion_ref;
          }

          // loop over z bins
          // copy content from ref phi bin, scaled
          for( int iz = 0; iz < h->GetNbinsZ(); ++iz )
          { h_copy->SetBinContent( ip+1, ir+1, iz+1, scale*h->GetBinContent( ip_ref+1, ir+1, iz+1 ) ); }
        }
      }
    }

    // write new histogram
    out->cd();
    h_copy->Write( hname );

  }
  out->Close();
}
