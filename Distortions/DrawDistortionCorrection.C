#include <RootUtil/Draw.h> 
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;
static constexpr double z_rec = 15;

// needed to determine z integration window
static constexpr double r_ref = 82;
static constexpr double z_ref = 33.25;
static constexpr double length = 50 - 5;
static constexpr double z_ref_min = z_ref - length/2;
static constexpr double z_ref_max = z_ref + length/2;

// z extrapolation between micromegas
// static constexpr double zextrap_min = 51.25;
// static constexpr double zextrap_max = 53.75;
static constexpr double zextrap_min = 48;
static constexpr double zextrap_max = 56;

namespace 
{ 
  int GetRangeBins( TAxis* axis ) { return axis->GetLast() - axis->GetFirst() + 1; }
}

//_______________________________________________
TString DrawDistortionCorrection()
{

  set_style( false );

  // open TFile  
  const TString tag = "_full_realistic_micromegas_mm-coarse";
  // const TString tag = "_full_realistic_micromegas_mm-coarse-new";

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

//   auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec_posz"));
//   auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec_posz"));
//   auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec_posz"));
//   auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec_posz"));

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
  std::cout << "DrawDistortionCorrection - z_rec: " << z_rec << std::endl;

  int phibin_rec = 0;
  double phi_rec_min = phi_rec;
  double phi_rec_max = phi_rec;

  int zbin_rec = 0;
  double z_rec_min = z_rec;
  double z_rec_max = z_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
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
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  int phibin_in = 0;
  int zbin_in_min = 0;
  int zbin_in_max = 0;
  int zbin_in = 0;
  
  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );
      phibin_in = h->GetXaxis()->FindBin( phi_rec );
      
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
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionCorrection - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in << ", " << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionCorrection - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in << ", " << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionCorrection%s_phi%i.pdf", tag.Data(), isec_rec ) );
  PdfDocument pdfDocument( pdfFile );

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

  if( hentries )
  {
    // entries in r, z plane
    TCanvas* cv( new TCanvas( "cv2", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    hentries->GetXaxis()->SetRange( phibin_rec, phibin_rec );
    hentries->GetZaxis()->SetRange( 0, 0 ); // z axis
    auto proj = hentries->Project3D( "yz" );
    proj->SetTitle("");
    proj->GetZaxis()->SetTitle( "entries" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw( "colz" );

    for( const auto& z:{z_ref, z_ref_min, z_ref_max} )
    {
      auto f = new TF1( "f", [](double*x, double *p){ return p[0] + p[1]*x[0]; }, 0, 100, 2 );
      f->SetParameter(0, 0);
      f->SetParameter(1, r_ref/z);
      f->SetLineColor(1);
      f->SetLineStyle(2);
      f->Draw( "same" );
    }
    
//     for( const auto& z:{zextrap_min, zextrap_max} )
//     {
//       auto f = new TF1( "f", [](double*x, double *p){ return p[0] + p[1]*x[0]; }, 0, 100, 2 );
//       f->SetParameter(0, 0);
//       f->SetParameter(1, r_ref/z);
//       f->SetLineColor(4);
//       f->SetLineStyle(2);
//       f->Draw( "same" );
//     }
    
    pdfDocument.Add(cv);

  }

  if( true )
  {
    // rdphi
    auto cv( new TCanvas( "cv3", "cv1", 800, 800 ) );
  file:///media/linux/hpereira/sphenix/work/g4simulations/macros/Distortions/DrawDistortionCorrection_z.C
  cv->SetLeftMargin( 0.15 );
    if( hDPint )
    {
      // deltaR vs R and z
      if( useRange ) 
      {
        hDPint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDPint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
      } else {
        hDPint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDPint->GetZaxis()->SetRange( zbin_in, zbin_in );
      }
      
      auto proj = hDPint->Project3D( "y" );
      proj->Scale( 1./(GetRangeBins(hDPint->GetXaxis())*GetRangeBins(hDPint->GetZaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "r (cm)" );
      proj->GetYaxis()->SetTitle( "r#Delta#phi (cm)" );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);
      
      if( addLimits )
      { 
        proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
        proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));      
      }
      
      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.65, Form( "#phi = %.3f rad", phi_rec ) );
      
//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionP_rec )
    {
      hDistortionP_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionP_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionP_rec->Project3D( "y" );
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
      
      if( useRange )
      {
        hDRint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDRint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
      } else {
        hDRint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDRint->GetZaxis()->SetRange( zbin_in, zbin_in );
      }
      
      // deltaR vs R and z
      auto proj = hDRint->Project3D( "y" );
      proj->Scale( 1./(GetRangeBins(hDRint->GetXaxis())*GetRangeBins(hDRint->GetZaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "r (cm)" );
      proj->GetYaxis()->SetTitle( "#Deltar (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      if( addLimits )
      { 
        proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
        proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));      
      }
      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.65, Form( "#phi = %.3f rad", phi_rec ) );

//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionR_rec )
    {
      hDistortionR_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionR_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionR_rec->Project3D( "y" );
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
    TCanvas* cv( new TCanvas( "cv5", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );

    if( hDZint )
    {
      if( useRange )
      {
        hDZint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
        hDZint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
      } else {
        hDZint->GetXaxis()->SetRange( phibin_in, phibin_in );
        hDZint->GetZaxis()->SetRange( zbin_in, zbin_in );        
      }
      
      // deltaR vs R and z
      auto proj = hDZint->Project3D( "y" );
      proj->Scale( 1./(GetRangeBins(hDZint->GetXaxis())*GetRangeBins(hDZint->GetZaxis())) );
      proj->SetTitle("");
      proj->GetXaxis()->SetTitle( "r (cm)" );
      proj->GetYaxis()->SetTitle( "#Deltaz (cm)" );
      proj->GetYaxis()->SetTitleOffset( 1.5 );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      // proj->SetMaximum( 0.4 );
      proj->SetMarkerStyle( 20 );
      proj->SetMarkerColor( 1 );
      proj->SetLineColor(1);

      if( addLimits )
      { 
        // proj->SetMaximum( 0.04 );
        proj->SetMinimum( std::min( proj->GetMinimum(), -0.1 ));
        proj->SetMaximum( std::max( proj->GetMaximum(), 0.1 ));      
      }
      proj->Draw( "P" );

      Draw::PutText( 0.6, 0.6, Form( "z = %.2f cm", z_rec ) );
      Draw::PutText( 0.6, 0.65, Form( "#phi = %.3f rad", phi_rec ) );

//       Draw::HorizontalLine( cv, 0 )->Draw();
//       Draw::VerticalLine( cv, 30 )->Draw();

    }

    if( hDistortionZ_rec )
    {
      hDistortionZ_rec->GetXaxis()->SetRange( phibin_rec, phibin_rec ); // phi axis
      hDistortionZ_rec->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
      auto proj = hDistortionZ_rec->Project3D( "y" );
      proj->SetMarkerStyle(20);
      proj->SetMarkerColor( 2 );
      proj->SetLineColor(2);
      proj->Draw( "same" );
    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
