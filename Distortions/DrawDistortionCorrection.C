#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr double z_rec = 10;

namespace
{
  int GetRangeBins( TAxis* axis )
  { return axis->GetLast() - axis->GetFirst() + 1; }

  constexpr double get_sector_phi( int isec )
  { return isec*M_PI/6; }
}

//_______________________________________________
TString DrawDistortionCorrection()
{

  set_style( false );
  gStyle->SetMarkerSize(1.);

  // const TString tag = "_offline_full_flat_acts_truth_distorted_truth";
  // const TString tag = "_offline_full_flat_genfit_truth_notpc_distorted_mm";
  // const TString tag = "_offline_full_flat_acts_truth_notpc_distorted_mm";
  // const TString tag = "_full_flat_acts_truth_notpc_distorted-new_mm";
  const TString tag = "_full_flat_genfit_truth_notpc_distorted_mm";
  const TString inputFile = Form("Rootfiles/Distortions%s.root", tag.Data());
  // const TString inputFile = Form("Distortions%s.root", tag.Data());

  const auto mapfile = "distortion_maps/average_minus_static_distortion_converted.root";
  // const auto mapfile = "distortion_maps/empty.root";
  TString pdfFile( Form( "Figures/DistortionCorrection%s.pdf", tag.Data() ) );

  bool use_phi_as_radian_input = true;

  std::cout << "DrawDistortionCorrection - input file: " << inputFile << std::endl;
  std::cout << "DrawDistortionCorrection - mapfile: " << mapfile << std::endl;
  std::cout << "DrawDistortionCorrection - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  const bool use_range = false;
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
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP_posz"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR_posz"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ_posz"));
  #endif

  if( hDistortionP_rec ) hDistortionP_rec->SetName( "hDistortionP_rec" );
  if( hDistortionR_rec ) hDistortionR_rec->SetName( "hDistortionR_rec" );
  if( hDistortionZ_rec ) hDistortionZ_rec->SetName( "hDistortionZ_rec" );

  if( hDistortionP_rec ) Utils::PrintAxis( hDistortionP_rec );
  else if( hDistortionZ_rec ) Utils::PrintAxis( hDistortionZ_rec );

  // find relevant reconstructed bins
  std::cout << "DrawDistortionCorrection - z_rec: " << z_rec << std::endl;

  int zbin_rec = 0;
  double z_rec_min = z_rec;
  double z_rec_max = z_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      zbin_rec = h->GetZaxis()->FindBin( z_rec );
      z_rec_min = h->GetZaxis()->GetBinLowEdge( zbin_rec );
      z_rec_max = h->GetZaxis()->GetBinUpEdge( zbin_rec );
      break;
    }
  }

  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get(z_rec > 0 ? "hIntDistortionR_posz":"hIntDistortionR_negz"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get(z_rec > 0 ? "hIntDistortionP_posz":"hIntDistortionP_negz"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get(z_rec > 0 ? "hIntDistortionZ_posz":"hIntDistortionZ_negz"));
  Utils::PrintAxis( hDPint );

  // find relevant input bins
  int zbin_in_min = 0;
  int zbin_in_max = 0;
  int zbin_in = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
      zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      zbin_in = h->GetZaxis()->FindBin( z_rec );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in << ", " << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  if( hentries )
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    hentries->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
    auto proj = hentries->Project3D( "yx" );
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

  if( hentries )
  {
    // entries in r, z plane
    TCanvas* cv( new TCanvas( "cv2", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    const int isec = 9;
    const auto phi_rec = get_sector_phi(isec);
    const auto phibin_rec = hentries->GetXaxis()->FindBin( phi_rec );

    hentries->GetXaxis()->SetRange( phibin_rec, phibin_rec );
    hentries->GetZaxis()->SetRange( 0, 0 ); // z axis
    auto proj = hentries->Project3D( "yz" );
    proj->SetTitle("");
    proj->GetZaxis()->SetTitle( "entries" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw( "colz" );
    pdfDocument.Add(cv);

  }

  // histogram pairs, for comparison
  using histogram_struct_t = std::tuple<TString, TH3*, TH3*, double>;
  std::array<histogram_struct_t, 3> histogram_structs = {{
    { "r#Delta#phi (cm)", hDPint, hDistortionP_rec, 0.5 },
    { "#Deltar (cm)", hDRint, hDistortionR_rec, 0.2 },
    { "#Deltaz (cm)", hDZint, hDistortionZ_rec, 0.04 }
  }};

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
      << "DrawDistortionCorrection -"
      << " isec: " << isec
      << " phibin_rec: " << phibin_rec
      << " zbin_rec: " << zbin_rec
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

    int ih =0;
    int ipad = 0;
    for( const auto& [label, hint, hDistortion_rec, limit]:histogram_structs )
    {

      cv->cd( ++ipad );
      gPad->SetLeftMargin( 0.15 );
      if( hint )
      {

        if( use_range )
        {
          hint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
          hint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
        } else {
          hint->GetXaxis()->SetRange( phibin_in, phibin_in );
          hint->GetZaxis()->SetRange( zbin_in, zbin_in );
        }

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
        hDistortion_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
        auto proj = hDistortion_rec->Project3D( "y" );

        // scale
        if( use_phi_as_radian_input && ih == 0 )
        {
          std::cout << "rescaling" << std::endl;
          for( int ix = 0; ix < proj->GetNbinsX(); ++ix )
          {
            const double r = proj->GetXaxis()->GetBinCenter( ix+1 );
            proj->SetBinContent( ix+1, proj->GetBinContent( ix+1 )*r );
          }
        }

        proj->SetMarkerStyle(20);
        proj->SetMarkerColor( 2 );
        proj->SetLineColor(2);
        proj->Draw( "same" );
      }

      ++ih;

    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
