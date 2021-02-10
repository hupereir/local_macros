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
static constexpr double r_rec = 70;

namespace 
{ 
  int GetRangeBins( TAxis* axis ) { return axis->GetLast() - axis->GetFirst() + 1; }
}

//_______________________________________________
TString DrawDistortionCorrection_z()
{

  /* bin definitions */
  set_style( false );

  // open TFile
  // const TString tag = "_full_realistic_micromegas_truth-coarse";
  const TString tag = "_full_realistic_micromegas_mm-coarse";
  // const TString tag = "_full_realistic_micromegas_mm-coarse-subtracted";
  // const TString tag = "_full_realistic_micromegas_mm-empty";
  // const TString tag = "_full_realistic_micromegas_truth-empty";
  const auto inputFile = Form( "Rootfiles/Distortions%s.root", tag.Data() );

  const bool addLimits = false;
  const bool useRange = false;
  
  auto f = TFile::Open( inputFile );
  std::cout << "DrawDistortionCorrection - input file: " << inputFile << std::endl;
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
  std::cout << "DrawDistortionCorrection - phi_rec: " << phi_rec << std::endl;
  std::cout << "DrawDistortionCorrection - r_rec: " << r_rec << std::endl;

  int phibin_rec = 0;
  double phi_rec_min = phi_rec;
  double phi_rec_max = phi_rec;

  int rbin_rec = 0;
  double r_rec_min = r_rec;
  double r_rec_max = r_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      phibin_rec = h->GetXaxis()->FindBin( phi_rec );
      phi_rec_min = h->GetXaxis()->GetBinLowEdge( phibin_rec );
      phi_rec_max = h->GetXaxis()->GetBinUpEdge( phibin_rec );

      rbin_rec = h->GetYaxis()->FindBin( r_rec );
      r_rec_min = h->GetYaxis()->GetBinLowEdge( rbin_rec );
      r_rec_max = h->GetYaxis()->GetBinUpEdge( rbin_rec );
      
      std::cout << "DrawDistortionCorrection - r_rec_min: " << r_rec_min << " max: " << r_rec_max << std::endl;

      break;
    }
  }

  // const auto mapfile = "distortion_maps/empty.root";
  const auto mapfile = "distortion_maps/fluct_average-coarse.root";
  // const auto mapfile = "distortion_maps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  // const auto mapfile = "distortion_maps/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ"));
  Utils::PrintAxis( hDPint );

  // find relevant input bins
  int phibin_in = 0;
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  int rbin_in = 0;
  int rbin_in_min = 0;
  int rbin_in_max = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );
      phibin_in = h->GetXaxis()->FindBin( phi_rec );
      
      rbin_in_min = h->GetYaxis()->FindBin( r_rec_min )+1;
      rbin_in_max = h->GetYaxis()->FindBin( r_rec_max );
      rbin_in = h->GetYaxis()->FindBin( r_rec );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in << ", " << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionCorrection - rbin_rec: " << rbin_rec << " rbin_in: (" << rbin_in << ", " << rbin_in_min << ", " << rbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionCorrection_z%s_phi%i.pdf", tag.Data(), isec_rec ) );
  PdfDocument pdfDocument( pdfFile );

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
      const double phi = (i+1)*M_PI/6;
      Draw::VerticalLine( cv, phi )->Draw();
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
    // rdphi
    auto cv( new TCanvas( "cv3", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );
    if( hDPint )
    {
      // deltaPhi vs z
      if( useRange )
      {
        hDPint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDPint->GetYaxis()->SetRange( rbin_in_min, rbin_in_max );
      } else {
        hDPint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDPint->GetYaxis()->SetRange( rbin_in, rbin_in );
      }
      
      auto proj = hDPint->Project3D( "z" );
      proj->Scale( 1./(GetRangeBins(hDPint->GetXaxis())*GetRangeBins(hDPint->GetYaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "z (cm)" );
      proj->GetYaxis()->SetTitle( "r#Delta#phi (cm)" );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
      proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));
      proj->Draw( "P" );

      Draw::PutText( 0.2, 0.8, Form( "r = %.2f cm", r_rec ) );
      Draw::PutText( 0.2, 0.85, Form( "#phi = %.3f rad", phi_rec ) );

//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionP_rec )
    {
      hDistortionP_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionP_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // r axis
      auto proj = hDistortionP_rec->Project3D( "z" );
      proj->SetMarkerStyle(20);
      proj->SetMarkerColor( 2 );
      proj->SetLineColor(2);
      proj->Draw( "same" );
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
      // deltaR vs z
      if( useRange )
      {
        hDRint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDRint->GetYaxis()->SetRange( rbin_in_min, rbin_in_max );
      } else {
        hDRint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDRint->GetYaxis()->SetRange( rbin_in, rbin_in );
      }
      
      auto proj = hDRint->Project3D( "z" );
      proj->Scale( 1./(GetRangeBins(hDRint->GetXaxis())*GetRangeBins(hDRint->GetYaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "z (cm)" );
      proj->GetYaxis()->SetTitle( "#Deltar (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
      proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));
      proj->Draw( "P" );

      Draw::PutText( 0.2, 0.8, Form( "r = %.2f cm", r_rec ) );
      Draw::PutText( 0.2, 0.85, Form( "#phi = %.3f rad", phi_rec ) );

//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionR_rec )
    {
      hDistortionR_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionR_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // z axis
      auto proj = hDistortionR_rec->Project3D( "z" );
      proj->SetMarkerStyle(20);
      proj->SetMarkerColor( 2 );
      proj->SetLineColor(2);
      proj->Draw( "same" );
    }

    pdfDocument.Add(cv);
  }

  if( true )
  {
    // dz corrections
    TCanvas* cv( new TCanvas( "cv5", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );

    if( hDZint )
    {
      if( useRange )
      {
        hDZint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDZint->GetYaxis()->SetRange( rbin_in_min, rbin_in_max );
      } else {
        hDZint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDZint->GetYaxis()->SetRange( rbin_in, rbin_in );
      }
      
      auto proj = hDZint->Project3D( "z" );
      proj->Scale( 1./(GetRangeBins(hDZint->GetXaxis())*GetRangeBins(hDZint->GetYaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "z (cm)" );
      proj->GetYaxis()->SetTitle( "#Deltaz (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      // proj->SetMaximum( 0.4 );
      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);
      proj->SetMinimum( -0.06 );

      proj->SetMinimum( std::min( proj->GetMinimum(), -0.04 ));
      proj->SetMaximum( std::max( proj->GetMaximum(), 0.04 ));
      proj->Draw( "P" );

      Draw::PutText( 0.2, 0.8, Form( "r = %.2f cm", r_rec ) );
      Draw::PutText( 0.2, 0.85, Form( "#phi = %.3f rad", phi_rec ) );

    }

    if( hDistortionZ_rec )
    {
      hDistortionZ_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionZ_rec->GetYaxis()->SetRange( rbin_rec, rbin_rec ); // z axis
      auto proj = hDistortionZ_rec->Project3D( "z" );
      proj->SetMarkerStyle(20);
      proj->SetMarkerColor( 2 );
      proj->SetLineColor(2);
      proj->Draw( "same" );
    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
