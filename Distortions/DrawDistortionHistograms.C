#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;
static constexpr double z_rec = 33.25;

// specify bins for which 2D histograms will be saved
static constexpr int phibin_rec = 11;
static constexpr int zbin_rec = 42;

//_______________________________________________
TString DrawDistortionHistograms()
{

  /* bin definitions */
  /**
  For Ross:
  PHG4TpcElectronDrift::InitRun - axis: phi bins: 38 limits: -0.174533 6.45772
  PHG4TpcElectronDrift::InitRun - axis: r bins: 18 limits: 16.375 81.625
  PHG4TpcElectronDrift::InitRun - axis: z bins: 42 limits: -2.6375 108.137

  For Me
  x axis: phi, 36 bins 0 to 2pi
  y axis: r, 16 bins, 20 to 78
  z axis: z, 80 bins -105.5, +105.5
  */

  set_style( false );

  // open TFile
  // const TString tag = "_full_realistic_micromegas_mm-old";
  const TString tag = "_full_realistic_micromegas_truth-empty";
  const auto inputfile = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  std::cout << "DrawDistortionHistograms - inputfile: " << inputfile << std::endl;
  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  // load histograms
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));

  // find relevant reconstructed bins
  std::cout << "DrawDistortionHistograms - phi_rec: " << phi_rec << std::endl;
  std::cout << "DrawDistortionHistograms - z_rec: " << z_rec << std::endl;

  int phibin_rec = 0;
  double phi_rec_min = phi_rec;
  double phi_rec_max = phi_rec;

  int zbin_rec = 0;
  double z_rec_min = z_rec;
  double z_rec_max = z_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      phibin_rec = h->GetXaxis()->FindBin( phi_rec );
      phi_rec_min = h->GetXaxis()->GetBinLowEdge( phibin_rec );
      phi_rec_max = h->GetXaxis()->GetBinUpEdge( phibin_rec );

      zbin_rec = h->GetZaxis()->FindBin( z_rec );
      z_rec_min = h->GetZaxis()->GetBinLowEdge( zbin_rec );
      z_rec_max = h->GetZaxis()->GetBinUpEdge( zbin_rec );
      break;
    }
  }

  const auto mapfile = "distortion_maps/empty.root";
  // const auto mapfile = "distortion_maps/fluct_average-coarse.root";
  // const auto mapfile = "distortion_maps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  // const auto mapfile = "distortion_maps/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ"));

  // find relevant input bins
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  int zbin_in_min = 0;
  int zbin_in_max = 0;

  for( const auto h: {hDPint, hDZint} )
  {
    if( h )
    {
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );

      zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
      zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionHistograms - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionHistograms - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionHistograms%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  // load phi histograms
  if( hDistortionP_rec )
  {

    auto cv = new TCanvas( "cv", "cv", 900, 900 );
    Draw::DivideCanvas( cv, hDistortionP_rec->GetYaxis()->GetNbins() );

    // project 3D dphi input histogram
    hDPint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDPint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDP = hDPint->Project3D( "y" );
    projDP->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // loop over rbins
    for( int ir = 0; ir < hDistortionP_rec->GetYaxis()->GetNbins(); ++ir )
    {
      const auto hname = Form( "residual_drphi_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH1*>( f->Get( hname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw();

      const auto value = hDistortionP_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      auto line = Draw::VerticalLine( gPad, value );
      line->SetLineColor( 2 );
      line->Draw();

      // find relevant input r bins
      {
        const auto r_min = hDistortionP_rec->GetYaxis()->GetBinLowEdge(ir+1);
        const auto r_max = hDistortionP_rec->GetYaxis()->GetBinUpEdge(ir+1);
        const auto rbin_min = projDP->GetXaxis()->FindBin( r_min );
        const auto rbin_max = projDP->GetXaxis()->FindBin( r_max );

        float dphi = 0;
        for( int irin = rbin_min; irin <= rbin_max; ++irin )
        { dphi += projDP->GetBinContent( irin ); }

        dphi /= (rbin_max-rbin_min+1);

        auto line = Draw::VerticalLine( gPad, dphi );
        line->SetLineColor( 4 );
        line->Draw();
      }
    }

    pdfDocument.Add(cv);

  }

  // load phi histograms
  if( hDistortionZ_rec )
  {

    auto cv = new TCanvas( "cv_z", "cv_z", 900, 900 );
    Draw::DivideCanvas( cv, hDistortionZ_rec->GetYaxis()->GetNbins() );

    // project 3D dphi input histogram
    hDZint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDZint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDZ = hDZint->Project3D( "y" );
    projDZ->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // loop over rbins
    for( int ir = 0; ir < hDistortionZ_rec->GetYaxis()->GetNbins(); ++ir )
    {
      const auto hname = Form( "residual_dz_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH1*>( f->Get( hname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw();

      const auto value = hDistortionZ_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      auto line = Draw::VerticalLine( gPad, value );
      line->SetLineColor( 2 );
      line->Draw();

      // find relevant input r bins
      {
        const auto r_min = hDistortionZ_rec->GetYaxis()->GetBinLowEdge(ir+1);
        const auto r_max = hDistortionZ_rec->GetYaxis()->GetBinUpEdge(ir+1);
        const auto rbin_min = projDZ->GetXaxis()->FindBin( r_min );
        const auto rbin_max = projDZ->GetXaxis()->FindBin( r_max );

        float dz = 0;
        for( int irin = rbin_min; irin <= rbin_max; ++irin )
        { dz += projDZ->GetBinContent( irin ); }

        dz /= (rbin_max-rbin_min+1);

        auto line = Draw::VerticalLine( gPad, dz );
        line->SetLineColor( 4 );
        line->Draw();
      }
    }

    pdfDocument.Add(cv);

  }
  return pdfFile;
}
