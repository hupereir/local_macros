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

static constexpr double z_rec = 5;
static constexpr double r_rec = 70;

//_______________________________________________
TString DrawDistortionCorrection_phi()
{

  set_style( false );

  // open TFile  
  const TString tag = "_full_realistic_micromegas_truth";

  auto f = TFile::Open( Form( "Rootfiles/Distortions%s.root", tag.Data() ) );
  // auto f = TFile::Open( Form( "distortion_maps_rec/Distortions%s.root", tag.Data() ) );
  if( !f ) return TString();

  // load histograms
  #if true
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
  #else
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ"));
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

  auto f2 = TFile::Open( "distortion_maps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" );
  // auto f2 = TFile::Open( "distortion_maps/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ"));
  Utils::PrintAxis( hDPint );

  // find relevant input bins
  int zbin_in_min = 0;
  int zbin_in_max = 0;

  int rbin_in_min = 0;
  int rbin_in_max = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      if( z_rec < 0 )
      {      
        zbin_in_min = h->GetZaxis()->FindBin( -z_rec_max )+1;
        zbin_in_max = h->GetZaxis()->FindBin( -z_rec_min );
      } else {
        zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
        zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      }

      rbin_in_min = h->GetYaxis()->FindBin( r_rec_min )+1;
      rbin_in_max = h->GetYaxis()->FindBin( r_rec_max );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionCorrection - rbin_rec: " << rbin_rec << " rbin_in: (" << rbin_in_min << ", " << rbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionCorrection_phi%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  if( true )
  {
    // rdphi
    auto cv( new TCanvas( "cv3", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );
    if( hDPint )
    {
      // deltaR vs R and z
      hDPint->GetYaxis()->SetRange( (rbin_in_min+rbin_in_max)/2, (rbin_in_min+rbin_in_max)/2 );
      hDPint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
      auto proj = hDPint->Project3D( "x" );
      proj->Scale( 1./(zbin_in_max-zbin_in_min+1) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "#phi (rad)" );
      proj->GetYaxis()->SetTitle( "r#Delta#phi (cm)" );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

//       proj->SetMinimum( -0.35 );
//       proj->SetMaximum( -0.3 );

      proj->SetMinimum( (proj->GetMinimum() > 0 ? 0.8:1.2)*proj->GetMinimum() );
      proj->SetMaximum( (proj->GetMaximum() > 0 ? 1.2:0.8)*proj->GetMaximum() );

      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.8, Form( "r = %.2f cm", r_rec ) );
      
    }

    if( hDistortionP_rec )
    {
      hDistortionP_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // r axis
      hDistortionP_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionP_rec->Project3D( "x" );
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
      
    pdfDocument.Add(cv);
  }

  if( true )
  {
    // dr corrections
    TCanvas* cv( new TCanvas( "cv4", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );
    if( hDRint )
    {
      hDRint->GetYaxis()->SetRange( (rbin_in_min+rbin_in_max)/2, (rbin_in_min+rbin_in_max)/2 );
      hDRint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );

      // deltaR vs R and z
      auto proj = hDRint->Project3D( "x" );
      proj->Scale( 1./(zbin_in_max-zbin_in_min+1) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "#phi (rad)" );
      proj->GetYaxis()->SetTitle( "#Deltar (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      proj->SetMinimum( (proj->GetMinimum() > 0 ? 0.8:1.2)*proj->GetMinimum() );
      proj->SetMaximum( (proj->GetMaximum() > 0 ? 1.2:0.8)*proj->GetMaximum() );

      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.8, Form( "r = %.2f cm", r_rec ) );

//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionR_rec )
    {
      hDistortionR_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // z axis
      hDistortionR_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionR_rec->Project3D( "x" );
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

    pdfDocument.Add(cv);
  }

  if( true )
  {
    // dr corrections
    TCanvas* cv( new TCanvas( "cv5", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );

    if( hDZint )
    {
      hDZint->GetYaxis()->SetRange( (rbin_in_min+rbin_in_max)/2, (rbin_in_min+rbin_in_max)/2 );
      hDZint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );

      // deltaR vs R and z
      auto proj = hDZint->Project3D( "x" );
      proj->Scale( 1./(zbin_in_max-zbin_in_min+1) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "#phi (rad)" );
      proj->GetYaxis()->SetTitle( "#Deltaz (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      proj->SetMinimum( (proj->GetMinimum() > 0 ? 0.8:1.2)*proj->GetMinimum() );
      proj->SetMaximum( (proj->GetMaximum() > 0 ? 1.2:0.8)*proj->GetMaximum() );

      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.8, Form( "r = %.2f cm", r_rec ) );

    }

    if( hDistortionZ_rec )
    {
      hDistortionZ_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // z axis
      hDistortionZ_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionZ_rec->Project3D( "x" );
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

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
