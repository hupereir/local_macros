#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
TString TrackHits()
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  static constexpr int nFiles = 2;
  const std::array<TString, nFiles> file =
  {
//     "DST/dst_eval_2k_realistic_full_nominal_old.root",
//     "DST/dst_eval_2k_realistic_full_nominal_new.root"
    "DST/dst_eval_full_high_occupancy_old.root",
    "DST/dst_eval_full_high_occupancy_new.root"
  };

  const TString pdfFile = "Figures/TrackHits_compare.pdf";

  static constexpr std::array<int, nFiles> color = { 1, 2 };
  static constexpr std::array<int, nFiles> symbol = { 20, 20 };
  const std::array<TString, nFiles> label =
  {
    "old errors",
    "new errors"
  };
 
  static constexpr int nvar = 4;
  const std::array<TString, nvar> var =
  {
    "_tracks._nclusters_mvtx", 
    "_tracks._nclusters_intt", 
    "_tracks._nclusters_tpc", 
    "_tracks._nclusters_ot",
  };

  const std::array<TString, nvar> varname =
  {
    "N_{clusters} MVTX", 
    "N_{clusters} INTT", 
    "N_{clusters} TPC", 
    "N_{clusters} OT" 
  };
  
  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  Draw::DivideCanvas( cv, nvar );

  // legends
  std::array<TLegend*, nvar> legend;
  for( int i = 0; i < nvar; ++i )
  { 
    legend[i] = new TLegend( 0.5, 0.8, 0.97, 0.9, "", "NDC" );
    legend[i]->SetFillColor(0);
    legend[i]->SetFillStyle(0);
    legend[i]->SetBorderSize(0);
  }
  
  // loop over files
  for( int i = 0; i < nFiles; ++i )
  {
    std::cout << "TrackHits - processing " << file[i] << std::endl;
    std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
    auto tree = static_cast<TTree*>( input->Get( "T" ) );
    gROOT->cd();
    
    // loop over variables
    for( int ivar = 0; ivar < nvar; ++ivar )
    {
      const auto hname = Form( "h%i%i", i, ivar );
      auto h = Utils::TreeToHisto( tree, hname, var[ivar], TCut(), true );
      h->SetTitle( "");
      h->GetXaxis()->SetTitle( varname[ivar] );
      h->SetLineColor( color[i] );
      cv->cd( ivar+1 );
      if( i ) h->Draw( "same" );
      else h->Draw();
      
      legend[ivar]->AddEntry( h, Form( "%s, mean=%.3f#pm%.3f", label[i].Data(), h->GetMean(), h->GetMeanError() ), "PL" );
      if( i == 0 ) 
      {
        legend[ivar]->Draw(); 
        gPad->SetLogy( true );
      }
    }
    
  }
  
  cv->SaveAs( pdfFile );
  return pdfFile;

}
