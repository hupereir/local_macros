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
void DrawClusterSize_Micromegas()
{
  
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 2;
  const std::array<TString,nFiles> file =
  {
    "Rootfiles/ClusterSize_Micromegas_Hijing_Micromegas_merged.root",
    "Rootfiles/ClusterSize_Micromegas_realistic_full_micromegas_nominal.root"
  };
  const std::array<TString,nFiles> labels = { "HIJING + PU", "Single Particles" };
  const TString pdfFile = "Figures/ClusterSize_Micromegas_Compare_Micromegas.pdf";

  constexpr std::array<int, nFiles> color = { 2, 1 };

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
    auto h55 = static_cast<TH1*>( input->Get( "profile_55" ) );
    auto h56 = static_cast<TH1*>( input->Get( "profile_56" ) );
    for( auto&& h:{h55,h56} )
    { 
      h->SetMaximum( 12 );
      h->SetLineColor( color[i] );
      h->SetMarkerColor( color[i] );
      h->SetMarkerStyle( 20 );
    }

    legend1->AddEntry( h55, labels[i], "LP" );
    legend2->AddEntry( h56, labels[i], "LP" );
    
    cv->cd(1);
    if( i ) h55->DrawCopy("same"); 
    else { h55->DrawCopy(); legend1->Draw(); }
    
    cv->cd(2);
    if( i ) h56->DrawCopy("same");
    else { h56->DrawCopy(); legend2->Draw(); }
  }
  
  cv->SaveAs( pdfFile );

}
