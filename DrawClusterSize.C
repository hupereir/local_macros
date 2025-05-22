#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
void DrawClusterSize()
{
  
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 3;
  const std::array<TString,nFiles> file =
  {
    "Rootfiles/ClusterSize_Hijing_micromegas_merged.root",
    "Rootfiles/ClusterSize_Hijing_micromegas_single.root",
    "Rootfiles/ClusterSize_1k_realistic_micromegas.root"
  };
  // const std::array<TString,nFiles> labels = { "HIJING + PU", "Single Particles" };
  const std::array<TString,nFiles> labels = { "HIJING + PU", "Hijing", "Single Particles" };
  const TString pdfFile = "Figures/ClusterSize_Compare_Micromegas.pdf";

  constexpr std::array<int, nFiles> color = { 2, 4, 1 };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 400 );
  cv->Divide( 2, 1 );
  
  auto legend1 = new TLegend( 0.5, 0.8, 0.97, 0.9, "", "NDC" );
  legend1->SetFillColor(0);
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);

  auto legend2 = new TLegend( 0.5, 0.8, 0.97, 0.9, "", "NDC" );
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(0);

  for( int i = 0; i < nFiles; ++i )
  {
    auto input( TFile::Open( file[i] ) );
    auto h55 = static_cast<TH1*>( input->Get( "h_55" ) );
    auto h56 = static_cast<TH1*>( input->Get( "h_56" ) );
    for( auto&& h:{h55,h56} )
    { 
      h->Scale( 1./h->GetEntries() );
      h->SetMaximum( 1 );
      h->SetLineColor( color[i] );
    }

    legend1->AddEntry( h55, Form( "%s, mean: %.1f", labels[i].Data(), h55->GetMean() ) , "L" );
    legend2->AddEntry( h56, Form( "%s, mean: %.1f", labels[i].Data(), h56->GetMean() ) , "L" );
    
    cv->cd(1);
    if( i ) h55->DrawCopy("hist same"); 
    else { h55->DrawCopy( "hist" ); legend1->Draw(); }
    gPad->SetLogy(true);
    
    cv->cd(2);
    if( i ) h56->DrawCopy("hist same");
    else { h56->DrawCopy( "hist" ); legend2->Draw(); }
    gPad->SetLogy(true);
  }
  
  cv->SaveAs( pdfFile );

}
