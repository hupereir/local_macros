#include <RootUtil/PdfDocument.h>
#include <RootUtil/Draw.h>

#include <TFile.h>
#include <TH3.h>

#include <array>
#include <cmath>

R__LOAD_LIBRARY(libRootUtilBase.so)

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

TString DistortionMapClosureTest()
{
  gStyle->SetOptStat(111111);

  const TString distortionfile = "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_converted.root";
  const TString correctionfile = "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_converted.root";
//  const TString correctionfile = "/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_inverted_10.root";
  const TString pdffile = "Figures/DistortionMapClosureTest.pdf";

  std::cout << "DistortionMapClosureTest - distortionfile: " << distortionfile << std::endl;
  std::cout << "DistortionMapClosureTest - correctionfile: " << correctionfile << std::endl;
  std::cout << "DistortionMapClosureTest - pdffile: " << pdffile << std::endl;

  PdfDocument pdfDocument( pdffile );

  // residual histograms
  const double dp_max = 0.001;
  auto hdp_diff = new TH1F( "hdp_diff", "", 100, -dp_max, dp_max );
  hdp_diff->GetXaxis()->SetTitle("#Delta#phi (rad)" );

  const double dr_max = 0.001;
  auto hdr_diff = new TH1F( "hdr_diff", "", 100, -dr_max, dr_max );
  hdr_diff->GetXaxis()->SetTitle("#Delta r (cm)" );

  const double dz_max = 0.001;
  auto hdz_diff = new TH1F( "hdz_diff", "", 100, -dz_max, dz_max );
  hdz_diff->GetXaxis()->SetTitle("#Delta z (cm)" );

  auto distortioninput = TFile::Open( distortionfile );
  auto correctioninput = TFile::Open( correctionfile );

  // loop over sides
  std::array<TString,2> suffix = { "negz", "posz" };
  for( int iside = 0; iside < 2; ++iside )
  {

    // load relevant histograms from distortioninput
    auto hdp_distortion = static_cast<TH3*>( distortioninput->Get( Form( "hIntDistortionP_%s", suffix[iside].Data() ) ) );
    auto hdr_distortion = static_cast<TH3*>( distortioninput->Get( Form( "hIntDistortionR_%s", suffix[iside].Data() ) ) );
    auto hdz_distortion = static_cast<TH3*>( distortioninput->Get( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) ) );
    print_histogram( hdp_distortion );

    // load relevant histograms from correctioninput
    auto hdp_correction = static_cast<TH3*>( correctioninput->Get( Form( "hIntDistortionP_%s", suffix[iside].Data() ) ) );
    auto hdr_correction = static_cast<TH3*>( correctioninput->Get( Form( "hIntDistortionR_%s", suffix[iside].Data() ) ) );
    auto hdz_correction = static_cast<TH3*>( correctioninput->Get( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) ) );
    print_histogram( hdp_correction );

    // loop over bins
    for( int ip = 0; ip < hdp_distortion->GetNbinsX()-2; ++ip )
      for( int ir = 0; ir < hdp_distortion->GetNbinsY()-2; ++ir )
      for( int iz = 0; iz < hdp_distortion->GetNbinsZ()-2; ++iz )
    {
      const double phi = hdp_distortion->GetXaxis()->GetBinCenter(ip+2);
      const double r = hdp_distortion->GetYaxis()->GetBinCenter(ir+2);
      const double z = hdp_distortion->GetZaxis()->GetBinCenter(iz+2);

      // get corresponding distortions
      const double dp_distortion = hdp_distortion->Interpolate(phi,r,z)/r;
      const double dr_distortion = hdr_distortion->Interpolate(phi,r,z);
      const double dz_distortion = hdz_distortion->Interpolate(phi,r,z);

      // get corrected position
      const double phi_distorted = phi+dp_distortion;
      const double r_distorted = r+dr_distortion;
      const double z_distorted = z+dz_distortion;

      // get corresponding correction
      const double dp_correction = hdp_correction->Interpolate(phi_distorted,r_distorted,z_distorted)/r_distorted;
      const double dr_correction = hdr_correction->Interpolate(phi_distorted,r_distorted,z_distorted);
      const double dz_correction = hdz_correction->Interpolate(phi_distorted,r_distorted,z_distorted);

      // apply correction
      const double phi_new = phi_distorted - dp_correction;
      const double r_new = r_distorted - dr_correction;
      const double z_new = z_distorted - dz_correction;

      // compare
      hdp_diff->Fill( phi_new - phi );
      hdr_diff->Fill( r_new - r );
      hdz_diff->Fill( z_new - z );
    }
  }


  // create plots
  for( const auto& h:{hdp_diff, hdr_diff, hdz_diff} )
  {
    const auto cvname = Form( "cv_%s", h->GetName() );
    auto cv = new TCanvas( cvname, cvname, 900, 900 );
    cv->SetRightMargin(0.15);
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );
    h->GetXaxis()->SetMaxDigits( 2 );
    h->GetYaxis()->SetMaxDigits( 4 );
    h->Draw("h");
    cv->Update();
    Draw::VerticalLine( cv, 0 )->Draw();

    pdfDocument.Add(cv);

  }

  return pdffile;

}
