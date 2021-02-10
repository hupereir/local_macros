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

TString DrawDistortionMap()
{

  set_style( false );

  // open TFile
  const TString tag = "_diff";
  const TString inputfile = "distortion_maps/diff_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";

  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  TString pdfFile( Form( "Figures/DistortionMap%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  // delta r vs r and z
  if( true )
  {

    // deltaR vs R and z
    auto h3= dynamic_cast<TH3*>(f->Get("hIntDistortionR"));
    h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
    Utils::PrintAxis( h3 );
    auto h = h3->Project3D( "yz" );
    h->SetTitle( "" );

    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    h->Scale( 1./h3->GetXaxis()->GetNbins() );
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "#Deltar (cm)" );
    h->GetZaxis()->SetTitleOffset( 1.7 );
    h->Draw("colz");
    gPad->SetRightMargin( 0.22 );
    Draw::HorizontalLine( cv, 30 )->Draw();
    pdfDocument.Add( cv );
  }

  // rdeltaphi vs r and z
  if( true )
  {

    // rdeltaphi vs R and z
    auto h3= dynamic_cast<TH3*>(f->Get("hIntDistortionP"));
    h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
    Utils::PrintAxis( h3 );
    auto h = h3->Project3D( "yz" );
    h->SetTitle( "" );

    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    h->Scale( 1./h3->GetXaxis()->GetNbins() );
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "r#Delta#phi (cm)" );
    h->GetZaxis()->SetTitleOffset( 1.7 );
    h->Draw("colz");
    gPad->SetRightMargin( 0.22 );
    Draw::HorizontalLine( cv, 30 )->Draw();
    pdfDocument.Add( cv );
  }

  // delta z vs r and z
  if( true )
  {

    // deltaR vs R and z
    auto h3= dynamic_cast<TH3*>(f->Get("hIntDistortionZ"));
    h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
    Utils::PrintAxis( h3 );
    auto h = h3->Project3D( "yz" );
    h->SetTitle( "" );

    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    h->Scale( 1./h3->GetXaxis()->GetNbins() );
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "#Deltaz (cm)" );
    h->GetZaxis()->SetTitleOffset( 1.7 );
    h->Draw("colz");
    gPad->SetRightMargin( 0.22 );
    Draw::HorizontalLine( cv, 30 )->Draw();
    pdfDocument.Add( cv );
  }

  return pdfFile;

}
