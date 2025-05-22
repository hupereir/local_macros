
#include <TFile.h>
#include <TH3.h>

#include <array>
#include <cmath>

namespace
{
  // print histogram
  void print_histogram( TH3* h )
  {

    std::cout << "InvertDistortions - name: " << h->GetName() << std::endl;
    for( const auto& axis:{h->GetXaxis(), h->GetYaxis(), h->GetZaxis() } )
    {
      std::cout
        << "  " << axis->GetName()
        << " bins: " << axis->GetNbins()
        << " min: " << axis->GetXmin()
        << " max: " << axis->GetXmax()
        << std::endl;
    }
    std::cout << std::endl;
  }

}

TString InvertDistortionMap()
{
//   const TString inputFile = "distortion_maps-new/average_minus_static_distortion_converted.root";
//   const TString outputFile = "distortion_maps-new/average_minus_static_distortion_inverted_10.root";
//   const bool use_phi_as_radian = false;
//   const int scale = 10;

  const TString inputFile = "/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX180961051_0.distortion_map.hist.root";
  const TString outputFile = "distortion_maps/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX180961051_0.distortion_map.inverted_10.root";
  const bool use_phi_as_radian = true;
  const int scale = 10;

  std::cout << "InvertDistortionMap - inputFile: " << inputFile << std::endl;
  std::cout << "InvertDistortionMap - outputFile: " << outputFile << std::endl;
  std::cout << "InvertDistortionMap - scale: " << scale << std::endl;

  auto input = TFile::Open( inputFile );
  auto output = TFile::Open( outputFile, "RECREATE" );

  std::array<TString,2> suffix = { "negz", "posz" };
  for( int iside = 0; iside < 2; ++iside )
  {

    // load relevan histograms from inputFile
    auto hdp = static_cast<TH3*>( input->Get( Form( "hIntDistortionP_%s", suffix[iside].Data() ) ) );
    auto hdr = static_cast<TH3*>( input->Get( Form( "hIntDistortionR_%s", suffix[iside].Data() ) ) );
    auto hdz = static_cast<TH3*>( input->Get( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) ) );

    print_histogram( hdp );

    // create converted histograms.
    /* they have different names to avoid conflicts */
    auto hdp_new = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionP_%s_inverted", suffix[iside].Data() ) ) );
    auto hdr_new = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionR_%s_inverted", suffix[iside].Data() ) ) );
    auto hdz_new = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionZ_%s_inverted", suffix[iside].Data() ) ) );

    // also create entries histogram
    auto h_entries = static_cast<TH3*>( hdp->Clone( Form( "entries_%s", suffix[iside].Data() ) ) );

    // reset
    for( const auto& h:{hdp_new, hdr_new, hdz_new, h_entries} ) { h->Reset(); }

    // get then number and limits of active bins
    const int npbins = hdp->GetNbinsX()-2;
    const auto p_min = hdp->GetXaxis()->GetBinLowEdge(2);
    const auto p_max = hdp->GetXaxis()->GetBinLowEdge(npbins+2);
    std::cout << "InvertDistortionMap - npbins: " << npbins << " p_min: " << p_min << " p_max: " << p_max << std::endl;

    const int nrbins = hdp->GetNbinsY()-2;
    const auto r_min = hdp->GetYaxis()->GetBinLowEdge(2);
    const auto r_max = hdp->GetYaxis()->GetBinLowEdge(nrbins+2);
    std::cout << "InvertDistortionMap - nrbins: " << nrbins << " r_min: " << r_min << " r_max: " << r_max << std::endl;

    const int nzbins = hdp->GetNbinsZ()-2;
    const auto z_min = hdp->GetZaxis()->GetBinLowEdge(2);
    const auto z_max = hdp->GetZaxis()->GetBinLowEdge(nzbins+2);
    std::cout << "InvertDistortionMap - nzbins: " << nzbins << " z_min: " << z_min << " z_max: " << z_max << std::endl;

    // scan inputs histograms
    for( int ip = 0; ip < npbins*scale; ++ip )
    {
      std::cout << "InvertDistortionMap ip: " << ip << "/" << npbins*scale << std::endl;
      for( int ir = 0; ir < nrbins*scale; ++ir )
      {
        for( int iz = 0; iz < nzbins*scale; ++iz )
        {

          // get initial position
          const double alpha_p = double(ip)/(npbins*scale-1);
          const double p = p_min*(1.-alpha_p) + p_max*alpha_p;

          const double alpha_r = double(ir)/(nrbins*scale-1);
          const double r = r_min*(1.-alpha_r) + r_max*alpha_r;

          const double alpha_z = double(iz)/(nzbins*scale-1);
          const double z = z_min*(1.-alpha_z) + z_max*alpha_z;

          // std::cout << "InvertDistortionMap - position: (" << p << ", " << r << ", " << z << ")" << std::endl;

          //  get distortions from source histograms
          const double dp = use_phi_as_radian ?
            (hdp->Interpolate(p,r,z)):
            (hdp->Interpolate(p,r,z)/r);
          const double dr = hdr->Interpolate(p,r,z);
          const double dz = hdz->Interpolate(p,r,z);

          // get new position
          double p_new = p+dp;
          if( p_new < 0 ) p_new += 2.*M_PI;
          if( p_new >= 2*M_PI ) p_new -= 2.*M_PI;

          const double r_new = r+dr;
          const double z_new = z+dz;

          // fill converted histogram
          if( use_phi_as_radian )
          {
            hdp_new->Fill( p_new, r_new, z_new, dp );
          } else {
            // important to multiply dp correction by r_new and not reuse drp
            const double drp_new = dp*r_new;
            hdp_new->Fill( p_new, r_new, z_new, drp_new );
          }

          hdr_new->Fill( p_new, r_new, z_new, dr );
          hdz_new->Fill( p_new, r_new, z_new, dz );
          h_entries->Fill( p_new, r_new, z_new );

        }
      }
    }

    // scale down by number of entries
    for( int ip = 0; ip < hdp_new->GetNbinsX(); ++ip )
      for( int ir = 0; ir < hdr_new->GetNbinsY(); ++ir )
      for( int iz = 0; iz < hdz_new->GetNbinsZ(); ++iz )
    {
      const double entries = h_entries->GetBinContent( ip+1, ir+1, iz+1 );
      std::cout << "InvertDistortionMap - (" << ip << ", " << ir << ", " << iz << ") - entries: " << entries << std::endl;
      if( entries > 0 )
      {
        for( const auto& h:{hdp_new, hdr_new, hdz_new} )
        { h->SetBinContent( ip+1, ir+1, iz+1, h->GetBinContent( ip+1, ir+1, iz+1 )/entries ); }
      }
    }

    // write to ouput TFile, using new name
    output->cd();
    hdp_new->Write( Form( "hIntDistortionP_%s", suffix[iside].Data() ) );
    hdr_new->Write( Form( "hIntDistortionR_%s", suffix[iside].Data() ) );
    hdz_new->Write( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) );

  }

  output->Close();
  return outputFile;
}
