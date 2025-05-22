#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr double z_rec = 1;

namespace
{
  int GetRangeBins( TAxis* axis ) 
  { return axis->GetLast() - axis->GetFirst() + 1; }
  
  constexpr double get_sector_phi( int isec ) 
  { return isec*M_PI/6; }
}

//_______________________________________________
/* this is the central membrane version for which distortion correction maps are 2D histograms instead of 3D */
TString DrawDistortionCorrection_cm()
{

  set_style( false );
  gStyle->SetMarkerSize(1.);
  
  // const TString tag = "-10ev";
  const TString tag = "_distorted-10ev-new";
  // const TString tag = "_distorted_scaled_x2-10ev-new";
  // const TString tag = "_distorted_full-10ev-new";
  const TString inputFile = Form("distortion_maps-new/CMDistortionCorrections%s.root", tag.Data());  

  const bool addLimits = true;

  std::cout << "DrawDistortionCorrection_cm - input file: " << inputFile << std::endl;
  auto f = TFile::Open( inputFile );
  if( !f ) return TString();

  // load histograms
  auto hentries = dynamic_cast<TH2*>(f->Get("hEntries_posz"));
  auto hDistortionP_rec = dynamic_cast<TH2*>(f->Get("hIntDistortionP_posz"));
  auto hDistortionR_rec = dynamic_cast<TH2*>(f->Get("hIntDistortionR_posz"));
  auto hDistortionZ_rec = dynamic_cast<TH2*>(f->Get("hIntDistortionZ_posz"));

  if( hDistortionP_rec ) hDistortionP_rec->SetName( "hDistortionP_rec" );
  if( hDistortionR_rec ) hDistortionR_rec->SetName( "hDistortionR_rec" );
  if( hDistortionZ_rec ) hDistortionZ_rec->SetName( "hDistortionZ_rec" );

  if( hDistortionP_rec ) Utils::PrintAxis( hDistortionP_rec );
  else if( hDistortionZ_rec ) Utils::PrintAxis( hDistortionZ_rec );

  // find relevant reconstructed bins
  std::cout << "DrawDistortionCorrection_cm - z_rec: " << z_rec << std::endl;
  
  // const auto mapfile = "distortion_maps/static_only.distortion_map.hist.root";
  const auto mapfile = "distortion_maps/average_minus_static_distortion_coarse.root";
  // const auto mapfile = "distortion_maps/average_minus_static_distortion_converted_scaled_x2.root";
  // const auto mapfile = "distortion_maps/empty.root";
  
  std::cout << "DrawDistortionCorrection_cm - mapfile: " << mapfile << std::endl;

  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR_posz"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP_posz"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ_posz"));
  Utils::PrintAxis( hDPint );

  // find relevant input bins
  int zbin_in = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      // always compare to positive z in the input map
      if( z_rec < 0 )
      {
        zbin_in = h->GetZaxis()->FindBin( -z_rec );
      } else {
        zbin_in = h->GetZaxis()->FindBin( z_rec );
      }
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection_cm - zbin_in: " << zbin_in <<std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionCorrection_cm%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  if( hentries )
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    auto proj = hentries;
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
      const auto phi_rec = get_sector_phi(isec);
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }

    pdfDocument.Add(cv);

  }

  // histogram pairs, for comparison
  using histogram_struct_t = std::tuple<TString, TH3*, TH2*, double>;
  std::array<histogram_struct_t, 3> histogram_structs = {{
    { "r#Delta#phi (cm)", hDPint, hDistortionP_rec, 0.5 },
    { "#Deltar (cm)", hDRint, hDistortionR_rec, 0.2 },
    { "#Deltaz (cm)", hDZint, hDistortionZ_rec, 0.04 }
  }};

//   using histogram_struct_t = std::tuple<TString, TH3*, TH2*, double>;
//   std::array<histogram_struct_t, 3> histogram_structs = {{
//     { "r#Delta#phi (cm)", hDPint, hDistortionP_rec, 0.7 },
//     { "#Deltar (cm)", hDRint, hDistortionR_rec, 0.3 },
//     { "#Deltaz (cm)", hDZint, hDistortionZ_rec, 0.06 }
//   }};

//     std::array<histogram_struct_t, 3> histogram_structs = {{
//     { "r#Delta#phi (cm)", hDPint, hDistortionP_rec, 2 },
//     { "#Deltar (cm)", hDRint, hDistortionR_rec, 2 },
//     { "#Deltaz (cm)", hDZint, hDistortionZ_rec, 0.1 }
//   }};


  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {
    const auto phi_rec = get_sector_phi(isec);
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

    std::cout 
      << "DrawDistortionCorrection_cm -"
      << " isec: " << isec 
      << " phibin_rec: " << phibin_rec 
      << std::endl;

    // find relevant input bins
    int phibin_in = 0;
    int phibin_in_min = 0;
    int phibin_in_max = 0;

    for( const auto h: {hDRint, hDPint, hDZint} )
    {
      if( h )
      {
        phibin_in = h->GetXaxis()->FindBin( phi_rec );
        phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
        phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );
        break;
      }
    }

    auto cv( new TCanvas( cvname, cvname, 900, 350 ) );
    cv->Divide( 3, 1 );

    int ipad = 0;
    for( const auto& [label, hint, hDistortion_rec, limit]:histogram_structs )
    {

      cv->cd( ++ipad );
      gPad->SetLeftMargin( 0.15 );
      if( hint )
      {
        hint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hint->GetZaxis()->SetRange( zbin_in, zbin_in );

        auto proj = hint->Project3D( "y" );
        proj->Scale( 1./(GetRangeBins(hint->GetXaxis())*GetRangeBins(hint->GetZaxis())) );
        proj->SetTitle("");
        proj->GetXaxis()->SetRangeUser( 30, 78 );
        proj->GetYaxis()->SetTitleOffset( 1.5 );
        proj->GetXaxis()->SetTitle( "r (cm)" );
        proj->GetYaxis()->SetTitle( label );

        // need to set errors to zero
        for( int i = 0; i < proj->GetNbinsX(); ++i )
        { proj->SetBinError( i+1, 0 ); }

        proj->SetMarkerStyle( 20 );
        proj->SetMarkerColor( 1 );
        proj->SetLineColor(1);

        if( addLimits )
        {
          proj->SetMinimum( std::min( proj->GetMinimum(), -limit));
          proj->SetMaximum( std::max( proj->GetMaximum(), limit ));
        }

        proj->Draw( "P" );

        Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
        Draw::PutText( 0.6, 0.65, Form( "#phi = %.3f rad", phi_rec ) );

      }

      if( hDistortion_rec )
      {
        hDistortion_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
        auto proj = hDistortion_rec->ProjectionY();
        proj->SetMarkerStyle(20);
        proj->SetMarkerColor( 2 );
        proj->SetLineColor(2);
        proj->Draw( "same" );
      }
    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
