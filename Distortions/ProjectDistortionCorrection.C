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

  // theta defined as atan2(z,r)
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

//______________________________________________________________________________
void ProjectDistortionCorrection()
{

  const TString inputfile = "distortion_maps/average_minus_static_distortion_converted.root";

  std::cout << "ProjectDistortionCorrection - inputfile: " << inputfile << std::endl;
  auto f = TFile::Open( inputfile );

  if( true )
  {
    // output filename
    const TString outputfile = "distortion_maps/average_minus_static_distortion_projected.root";
    std::cout << "ProjectDistortionCorrection - outputfile: " << outputfile << std::endl;
    auto out = TFile::Open( outputfile, "RECREATE" );

    auto hentries_posz = static_cast<TH3*>( f->Get("hIntDistortionR_posz")->Clone("hentries_posz"));
    auto hentries_negz = static_cast<TH3*>( f->Get("hIntDistortionR_negz")->Clone("hentries_negz"));
    hentries_posz->Reset();
    hentries_negz->Reset();

    for( const TString& hname: {
      "hIntDistortionR_posz",  "hIntDistortionR_negz",
      "hIntDistortionP_posz",  "hIntDistortionP_negz",
      "hIntDistortionZ_posz",  "hIntDistortionZ_negz" } )
    {
      // get histogram
      auto h = static_cast<TH3*>( f->Get(hname) );

      // clone
      auto h_copy = static_cast<TH3*>(h->Clone(hname+"_copy"));
      h_copy->Reset();

      // loop over bins
      for( int ip = 0; ip < h->GetNbinsX(); ++ip )
        for( int ir = 0; ir < h->GetNbinsY(); ++ir )
        for( int iz = 0; iz < h->GetNbinsZ(); ++iz )
      {
        const double phi = h->GetXaxis()->GetBinCenter(ip+1);
        const double r = h->GetYaxis()->GetBinCenter(ir+1);
        const double z = h->GetZaxis()->GetBinCenter(iz+1);
        const double theta = std::atan2( z, r );

        if(
          ( in_range( phi, phi_window_central ) && in_range( theta, theta_window_central ) ) ||
          ( in_range( phi, phi_window_east ) && in_range( theta, theta_window_east ) ) ||
          ( in_range( phi, phi_window_west ) && in_range( theta, theta_window_west ) ))
        {
          hentries_posz->Fill( phi, r, z );
          hentries_negz->Fill( phi, r, z );
          h_copy->SetBinContent( ip+1, ir+1, iz+1, h->GetBinContent(ip+1, ir+1, iz+1) );
        }
      }

      out->cd();
      h_copy->SetName( hname );
      h_copy->Write( hname );

    }

    out->cd();
    hentries_posz->Write();
    hentries_negz->Write();
    out->Close();
  }

  if( true )
  {
    // get central membrane 2D histograms
    const TString outputfile_cm = "distortion_maps/average_minus_static_distortion_cm.root";
    std::cout << "ProjectDistortionCorrection - outputfile_cm: " << outputfile_cm << std::endl;
    auto out_cm = TFile::Open( outputfile_cm, "RECREATE" );

    // get 2D histogram from 3D
    auto project_to_cm = []( TH3* h, int bin )
    {
      h->GetZaxis()->SetRange( bin, bin );
      return h->Project3D( "yx" );
    };

    // posz histograms
    for( const TString& hname: {
      "hIntDistortionR_posz",
      "hIntDistortionP_posz",
      "hIntDistortionZ_posz" } )
    {
      // get histogram
      auto h = static_cast<TH3*>( f->Get(hname) );

      // project the first bin (z==0) into 2D histogram
      auto h_2d = project_to_cm( h, 1 );
      out_cm->cd();
      h_2d->SetName( hname );
      h_2d->Write( hname );
    }

    // negz histograms
    for( const TString& hname: {
      "hIntDistortionR_negz",
      "hIntDistortionP_negz",
      "hIntDistortionZ_negz" } )
    {
      // get histogram
      auto h = static_cast<TH3*>( f->Get(hname) );

      // project the last bin (z==0) into 2D histogram
      auto h_2d = project_to_cm( h, h->GetNbinsZ() );
      out_cm->cd();
      h_2d->SetName( hname );
      h_2d->Write( hname );
    }
    out_cm->Close();

  }

}
