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

TString DrawDistortionMap_rphi()
{

  set_style( false );

  // open TFile
  // auto f = TFile::Open( "distortion_maps/BeamXingNBeamsx10.flat_B1.4_E-400.0.ross_phislice_lookup_r16xp36xz40.distortion_map.hist.root" );
  // auto f = TFile::Open( "distortion_maps/average.rev1.hist0.realE_B-1.5_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" );

  const TString tag = "averaged";
  auto f = TFile::Open( "distortion_maps/output.averaged.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" );

//   const TString tag = "single";
//   auto f = TFile::Open( "distortion_maps/output.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" );
  if( !f ) return TString();

  TString pdfFile = Form( "Figures/distortionMap_%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f->Get("hIntDistortionR"));
  auto hDPint= dynamic_cast<TH3*>(f->Get("hIntDistortionP"));

  if( true )
  {
    TCanvas* cv( new TCanvas( "cv0", "cv0", 800, 800 ) );
    cv->SetRightMargin( .15 );

    // z bin
    hDPint->GetZaxis()->SetRange( 2, 2 );
    auto proj = hDPint->Project3D( "yx" );
    proj->SetTitle( "" );
    proj->GetZaxis()->SetTitle( "r#Delta#phi" );
    proj->Draw( "colz" );
    pdfDocument.Add( cv );
  }

  if( true )
  {

    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );

    for( int iphi = 2; iphi < 38; ++iphi )
    {

      hDPint->GetXaxis()->SetRange( iphi, iphi );
      hDPint->GetZaxis()->SetRange( 2, 2 );

      // deltaR vs R and z
      auto proj = hDPint->Project3D( "y" );
      proj->SetName( Form( "proj_%i", iphi ) );
      proj->SetTitle( "" );

      // need to set errors to zero
      for( int i = 0; i < proj->GetNbinsX(); ++i )
      { proj->SetBinError( i+1, 0 ); }

      proj->SetMarkerStyle( 20 );
      proj->GetYaxis()->SetTitle( "r#Delta#phi" );
      if( iphi == 2 ) proj->Draw();
      else proj->Draw("same" );
    }

    pdfDocument.Add( cv );
  }

  return pdfFile;

}
