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

void DrawDeltaRPhi_pt()
{
  
  set_style( false );
  
  const TString tag =  "_flat_full_micromegas_notpc-nomicromegas-new" ;
  const TString inputfile( Form( "Rootfiles/DeltaRPhi_truth%s.root", tag.Data() ) );
  const TString pdffile( Form( "Figures/DeltaRPhi_truth_pt%s.pdf", tag.Data() ) );
  
  auto input = TFile::Open( inputfile );
   
  constexpr int nptbins = 20;
  std::array<float, nptbins+1> invptbins = {{ 
    0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 
    1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05 }};

  constexpr std::array<int, nptbins> color = {{
    kMagenta+1, kBlue+2, kBlue, kCyan+2, kGreen+2, kGreen+1, kOrange, kOrange-3, kRed+2, kRed,
    kMagenta+1, kBlue+2, kBlue, kCyan+2, kGreen+2, kGreen+1, kOrange, kOrange-3, kRed+2, kRed
  }};

  auto cv = new TCanvas( "cvtg", "cvtg", 800, 600 );
  cv->SetLeftMargin( 0.16 );

  auto h = new TH1F( "dummy", "", 100, 20, 80 );
  h->SetMinimum(0);
  h->SetMaximum( 8000 );
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.6 );
  h->Draw();
  
  auto legend = new TLegend( 0.15, 0.70, 0.50, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();
  
  auto legend2 = new TLegend( 0.45, 0.70, 0.80, 0.9, "", "NDC" );
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(0);
  legend2->Draw();
  
  for( int ibin = 0; ibin < 10; ++ibin )
  {
    
    TString tgname = Form( "residuals_%i", ibin );
    auto tg = static_cast<TGraphErrors*>(  input->Get( tgname ) );
    if( !tg ) continue;
    
    tg->SetMarkerStyle(20);
    tg->SetLineColor(color[ibin]);
    tg->SetMarkerColor(color[ibin]);

    if( ibin < 5 ) legend->AddEntry( tg, Form( "1/#it{p}_{T} = %.1f GeV^{-1}", 0.5*(invptbins[ibin]+invptbins[ibin+1]) ), "AP" );
    else legend2->AddEntry( tg, Form( "1/#it{p}_{T} = %.1f GeV^{-1}", 0.5*(invptbins[ibin]+invptbins[ibin+1]) ), "AP" );
    tg->Draw( "P" );
    
  }
  
  cv->SaveAs( pdffile );
  
}
