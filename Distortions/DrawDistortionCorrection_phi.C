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
static constexpr double phi_rec = isec_rec*M_PI/6 + M_PI/12;

static constexpr double z_rec = 15;
static constexpr double r_rec = 60;

namespace
{
  int GetRangeBins( TAxis* axis ) { return axis->GetLast() - axis->GetFirst() + 1; }
}

//_______________________________________________
TString DrawDistortionCorrection_phi()
{

  set_style( false );


  const TString tag = "_extrapolated";
  const TString inputFile = "distortion_maps/average_minus_static_distortion_extrapolated.root";
  const auto mapfile = "distortion_maps/average_minus_static_distortion_converted.root";
  // const auto mapfile = "distortion_maps/empty.root";

  const TString pdfFile( Form( "Figures/DistortionCorrection_phi%s.pdf", tag.Data() ) );

  std::cout << "DrawDistortionCorrection - input file: " << inputFile << std::endl;
  std::cout << "DrawDistortionCorrection - mapfile: " << mapfile << std::endl;
  std::cout << "DrawDistortionCorrection - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  const bool addLimits = false;
  const bool useRange = false;

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
  std::cout << "DrawDistortionCorrection - r_rec: " << r_rec << std::endl;
  std::cout << "DrawDistortionCorrection - z_rec: " << z_rec << std::endl;

  int zbin_rec = 0;
  double z_rec_min = z_rec;
  double z_rec_max = z_rec;

  int rbin_rec = 0;
  double r_rec_min = r_rec;
  double r_rec_max = r_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      zbin_rec = h->GetZaxis()->FindBin( z_rec );
      z_rec_min = h->GetZaxis()->GetBinLowEdge( zbin_rec );
      z_rec_max = h->GetZaxis()->GetBinUpEdge( zbin_rec );

      rbin_rec = h->GetYaxis()->FindBin( r_rec );
      r_rec_min = h->GetYaxis()->GetBinLowEdge( rbin_rec );
      r_rec_max = h->GetYaxis()->GetBinUpEdge( rbin_rec );
      break;
    }
  }

  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR_posz"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP_posz"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ_posz"));
  Utils::PrintAxis( hDPint );

  // find relevant input bins
  int zbin_in_min = 0;
  int zbin_in_max = 0;
  int zbin_in = 0;

  int rbin_in_min = 0;
  int rbin_in_max = 0;
  int rbin_in = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      if( z_rec < 0 )
      {
        zbin_in_min = h->GetZaxis()->FindBin( -z_rec_max )+1;
        zbin_in_max = h->GetZaxis()->FindBin( -z_rec_min );
        zbin_in = h->GetZaxis()->FindBin( -z_rec );
      } else {
        zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
        zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
        zbin_in = h->GetZaxis()->FindBin( z_rec );
      }

      rbin_in_min = h->GetYaxis()->FindBin( r_rec_min )+1;
      rbin_in_max = h->GetYaxis()->FindBin( r_rec_max );
      rbin_in =  h->GetYaxis()->FindBin( r_rec );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionCorrection - rbin_rec: " << rbin_rec << " rbin_in: (" << rbin_in_min << ", " << rbin_in_max << ")" << std::endl;

  // histogram pairs, for comparison
  using histogram_struct_t = std::tuple<TString, TH3*, TH3*>;
  std::array<histogram_struct_t, 3> histogram_pairs = {{
    { "r#Delta#phi (cm)", hDPint, hDistortionP_rec },
    { "#Deltar (cm)", hDRint, hDistortionR_rec },
    { "#Deltaz (cm)", hDZint, hDistortionZ_rec }
  }};

  auto cv( new TCanvas( "cv", "cv", 900, 350 ) );
  cv->Divide( 3, 1 );
  int cv_id = 0;

  for( const auto& [label, hint, hDistortion_rec]:histogram_pairs )
  {
    cv->cd(++cv_id);
    gPad->SetLeftMargin( 0.15 );
    if( hint )
    {
      if( useRange )
      {
        hint->GetYaxis()->SetRange( rbin_in_min, rbin_in_max );
        hint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
      } else {
        hint->GetYaxis()->SetRange( rbin_in, rbin_in );
        hint->GetZaxis()->SetRange( zbin_in, zbin_in );
      }

      auto proj = hint->Project3D( "x" );
      proj->Scale( 1./(GetRangeBins(hint->GetYaxis())*GetRangeBins(hint->GetZaxis())) );
      proj->SetTitle("");
      proj->GetYaxis()->SetTitleOffset( 1.5 );
      proj->GetYaxis()->SetTitle( label );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

//       if( addLimits )
//       {
//         proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
//         proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));
//       }

      proj->SetMinimum( (proj->GetMinimum() > 0 ? 0.8:1.2)*proj->GetMinimum() );
      proj->SetMaximum( (proj->GetMaximum() > 0 ? 1.2:0.8)*proj->GetMaximum() );

      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.8, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.85, Form( "r = %.2f cm", r_rec ) );

    }

    if( hDistortion_rec )
    {
      hDistortion_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // r axis
      hDistortion_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortion_rec->Project3D( "x" );
      proj->SetMarkerStyle(20);
      proj->SetMarkerColor( 2 );
      proj->SetLineColor(2);
      proj->Draw( "same" );
    }

    {
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }

  }

  pdfDocument.Add(cv);
  return pdfFile;
}
