#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libtpccalib.so)

// static constexpr double r_rec = 43.56;
// static constexpr double r_rec = 60;
// static constexpr double r_rec = 70;
static constexpr double r_rec = 40;

namespace
{
  int get_range_bins( TAxis* axis ) { return axis->GetLast() - axis->GetFirst() + 1; }
  double get_sector_phi( int isec ) { return isec*M_PI/6; }
}

//_______________________________________________
TString DrawDistortionCorrection_z()
{

  /* bin definitions */
  set_style( false );

  // open TFile
  // const TString tag = "_full_flat_acts_truth_no_distortion_truth";
  // const TString tag = "_full_flat_acts_truth_distorted_truth";
  // const TString tag = "_offline_full_flat_acts_truth_notpc_nodistortion_mm";
  // const TString tag = "_offline_full_flat_genfit_truth_notpc_nodistortion_mm";
  //   const TString tag = "_offline_full_flat_genfit_truth_notpc_nodistortion-nosurvey_mm";
  //   const TString inputFile = Form("Rootfiles/Distortions%s.root", tag.Data());

  const TString tag = "_extrapolated";
  const TString inputFile = "distortion_maps/average_minus_static_distortion_extrapolated.root";
  const auto mapfile = "distortion_maps/average_minus_static_distortion_converted.root";
  // const auto mapfile = "distortion_maps/empty.root";

  const TString pdfFile( Form( "Figures/DistortionCorrection_z%s.pdf", tag.Data() ) );

  std::cout << "DrawDistortionCorrection - input file: " << inputFile << std::endl;
  std::cout << "DrawDistortionCorrection - mapfile: " << mapfile << std::endl;
  std::cout << "DrawDistortionCorrection - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  const bool useRange = false;
  const bool addLimits = true;

  auto f = TFile::Open( inputFile );
  if( !f ) return TString();

  // load histograms
  #if false
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
  #else
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP_negz"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR_negz"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ_negz"));
//   auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP_posz"));
//   auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR_posz"));
//   auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ_posz"));
  #endif

  if( hDistortionP_rec ) hDistortionP_rec->SetName( "hDistortionP_rec" );
  if( hDistortionR_rec ) hDistortionR_rec->SetName( "hDistortionR_rec" );
  if( hDistortionZ_rec ) hDistortionZ_rec->SetName( "hDistortionZ_rec" );

  if( hDistortionP_rec ) Utils::PrintAxis( hDistortionP_rec );
  else if( hDistortionZ_rec ) Utils::PrintAxis( hDistortionZ_rec );

  // find relevant reconstructed bins
  std::cout << "DrawDistortionCorrection - r_rec: " << r_rec << std::endl;

  int rbin_rec = 0;
  double r_rec_min = r_rec;
  double r_rec_max = r_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      rbin_rec = h->GetYaxis()->FindBin( r_rec );
      r_rec_min = h->GetYaxis()->GetBinLowEdge( rbin_rec );
      r_rec_max = h->GetYaxis()->GetBinUpEdge( rbin_rec );
      std::cout << "DrawDistortionCorrection - r_rec_min: " << r_rec_min << " max: " << r_rec_max << std::endl;

      break;
    }
  }

  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  std::array<TH3*,2> hDRint = { dynamic_cast<TH3*>(f2->Get("hIntDistortionR_negz")), dynamic_cast<TH3*>(f2->Get("hIntDistortionR_posz"))};
  std::array<TH3*,2> hDPint = { dynamic_cast<TH3*>(f2->Get("hIntDistortionP_negz")), dynamic_cast<TH3*>(f2->Get("hIntDistortionP_posz"))};
  std::array<TH3*,2> hDZint = { dynamic_cast<TH3*>(f2->Get("hIntDistortionZ_negz")), dynamic_cast<TH3*>(f2->Get("hIntDistortionZ_posz"))};

  Utils::PrintAxis( hDPint[0] );

  // find relevant input bins
  int rbin_in = 0;
  int rbin_in_min = 0;
  int rbin_in_max = 0;

  for( const auto h: {hDRint[0], hDPint[0], hDZint[0]} )
  {
    if( h )
    {
      rbin_in_min = h->GetYaxis()->FindBin( r_rec_min )+1;
      rbin_in_max = h->GetYaxis()->FindBin( r_rec_max );
      rbin_in = h->GetYaxis()->FindBin( r_rec );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - rbin_rec: " << rbin_rec << " rbin_in: (" << rbin_in << ", " << rbin_in_min << ", " << rbin_in_max << ")" << std::endl;

  if( hentries )
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    hentries->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // r axis
    auto proj = hentries->Project3D( "zx" );
    proj->SetTitle("");
    proj->GetZaxis()->SetTitle( "entries" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw( "colz" );

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const auto phi( get_sector_phi(i) + M_PI/12 );
      Draw::VerticalLine( cv, phi )->Draw();
    }

    {
      const int isec = 9;
      const auto phi_rec = get_sector_phi( isec );
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }

    pdfDocument.Add(cv);

  }

  // histogram pairs, for comparison
  using histogram_struct_t = std::tuple<TString, TH3*, TH3*, double>;
//   std::array<histogram_struct_t, 3> histogram_pairs = {{
//     { "r#Delta#phi (cm)", hDPint[1], hDistortionP_rec, 0.1 },
//     { "#Deltar (cm)", hDRint[1], hDistortionR_rec, 0.1 },
//     { "#Deltaz (cm)", hDZint[1], hDistortionZ_rec, 0.1 }
//   }};

  std::array<histogram_struct_t, 3> histogram_pairs = {{
    { "r#Delta#phi (cm)", hDPint[0], hDistortionP_rec, 0.1 },
    { "#Deltar (cm)", hDRint[0], hDistortionR_rec, 0.1 },
    { "#Deltaz (cm)", hDZint[0], hDistortionZ_rec, 0.1 }
  }};

//   // get reference z range
//   const auto [zmin,zmax] = TpcSpaceChargeReconstructionHelper::get_zref_range( r_rec );

  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {
    const auto phi_rec = get_sector_phi( isec );
    const auto cvname = Form( "cv%i", isec );

    int phibin_rec = 0;
    double phi_rec_min = phi_rec;
    double phi_rec_max = phi_rec;

    for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
    {
      if( h )
      {
        phibin_rec = h->GetXaxis()->FindBin( phi_rec );
        phi_rec_min = h->GetXaxis()->GetBinLowEdge( phibin_rec );
        phi_rec_max = h->GetXaxis()->GetBinUpEdge( phibin_rec );
        break;
      }
    }

    // find relevant input bins
    int phibin_in = 0;
    int phibin_in_min = 0;
    int phibin_in_max = 0;

    for( const auto h: {hDRint[0], hDPint[0], hDZint[0]} )
    {
      if( h )
      {
        phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
        phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );
        phibin_in = h->GetXaxis()->FindBin( phi_rec );
        break;
      }
    }

    auto cv( new TCanvas( cvname, cvname, 900, 350 ) );
    cv->Divide( 3, 1 );

    int ipad = 0;
    for( const auto& [label, hint, hDistortion_rec, limit]:histogram_pairs )
    {

      cv->cd( ++ipad );
      gPad->SetLeftMargin( 0.15 );
      if( hint )
      {
        // deltaPhi vs z
        if( useRange )
        {
          hint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
          hint->GetYaxis()->SetRange( rbin_in_min, rbin_in_max );
        } else {
          hint->GetXaxis()->SetRange( phibin_in, phibin_in );
          hint->GetYaxis()->SetRange( rbin_in, rbin_in );
        }

        auto proj = hint->Project3D( "z" );
        proj->Scale( 1./(get_range_bins(hint->GetXaxis())*get_range_bins(hint->GetYaxis())) );
        proj->SetTitle("");
        proj->GetYaxis()->SetTitleOffset( 1.5 );
        proj->GetYaxis()->SetTitle( label );

        // need to set errors to zero
        for( int i = 0; i < proj->GetNbinsX(); ++i )
        { proj->SetBinError( i+1, 0 ); }

        proj->SetMarkerStyle( 20 );
        proj->SetMarkerColor( 1 );
        proj->SetLineColor(1);

        proj->SetMinimum( std::min( proj->GetMinimum(), -limit ));
        proj->SetMaximum( std::max( proj->GetMaximum(), limit ));
        proj->Draw( "P" );

        Draw::PutText( 0.2, 0.8, Form( "r = %.2f cm", r_rec ) );
        Draw::PutText( 0.2, 0.85, Form( "#phi = %.3f rad", phi_rec ) );

      }

      if( hDistortion_rec )
      {
        hDistortion_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
        hDistortion_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // r axis
        auto proj = hDistortion_rec->Project3D( "z" );
        proj->SetMarkerStyle(20);
        proj->SetMarkerColor( 2 );
        proj->SetLineColor(2);
        proj->Draw( "same" );

//         for( const auto z:{zmin, zmax} )
//         { Draw::VerticalLine( gPad, z )->Draw(); }


//         // get relevant integral
//         {
//           const auto zbinmin = hDistortion_rec->GetZaxis()->FindBin( zmin );
//           const auto zbinmax = hDistortion_rec->GetZaxis()->FindBin( zmax );
//           const auto norm = hDistortion_rec->Integral( phibin_rec, phibin_rec, rbin_rec, rbin_rec, zbinmin, zbinmax ) / (zbinmax-zbinmin+1);
//           TLine *line = new TLine( zmin, norm, zmax, norm );
//           line->SetLineStyle( 2 );
//           line->SetLineWidth( 2 );
//           line->SetLineColor( 4 );
//           line->Draw();
//         }

      }

    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
