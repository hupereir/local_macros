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
void DrawDeltaRPhi()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
// 
//   constexpr int nfiles = 2;
//   const std::array<TString, nfiles> file =
//   {
// //     "Rootfiles/DeltaRPhi_upsilon_acts_full_no_distortion-new.root",
// //     "Rootfiles/DeltaRPhi_upsilon_acts_full_distorted-new.root"
// //     "Rootfiles/DeltaRPhi_upsilon_acts_truth_no_distortion-new.root",
// //     "Rootfiles/DeltaRPhi_upsilon_acts_truth_distorted-new.root"
// 
//     "Rootfiles/DeltaRPhi_upsilon_acts_truth_no_distortion-new.root",
//     "Rootfiles/DeltaRPhi_upsilon_acts_truth_distorted_fullmap-new.root"
//   };
// 
//   // const TString pdfFile = "Figures/DeltaRPhi_compare_full_distortion-new.pdf";
//   // const TString pdfFile = "Figures/DeltaRPhi_compare_truth_distortion-new.pdf";
//   
//   constexpr std::array<int, 2> color = { kBlack, kRed };
//   constexpr std::array<int, 2> symbol = { 20, 20 };
//   const std::array<TString, 2> labels = 
//   { 
//     "w/o distortions",
//     // "w/ beam-induced distortions"
//     "w/ static distortions"
//   };
// 
//   const TString pdfFile = "Figures/DeltaRPhi_compare_truth_distortion_fullmap-new.pdf";

//   constexpr int nfiles = 3;
//   const std::array<TString, nfiles> file =
//   {
//     "Rootfiles/DeltaRPhi_upsilon_acts_full_no_distortion-new.root",
//     "Rootfiles/DeltaRPhi_upsilon_acts_full_distorted-tony.root",
//     "Rootfiles/DeltaRPhi_upsilon_acts_full_distorted_fullmap-tony.root"
//   };
// 
//   // const TString pdfFile = "Figures/DeltaRPhi_compare_full_distortion-new.pdf";
//   // const TString pdfFile = "Figures/DeltaRPhi_compare_truth_distortion-new.pdf";
//   
//   constexpr std::array<int, 3> color = { 1, 2, 4 };
//   constexpr std::array<int, 3> symbol = { 20, 20, 20 };
// 
//   std::array<TString, nfiles> labels =
//   {
//     "no distortions",
//     "beam-induced distortions + correction",
//     "static distortions + correction"
//   };
//   
//   const TString pdfFile = "Figures/DeltaRPhi_compare_full_distortions.pdf";

  
    constexpr int nfiles = 3;
  const std::array<TString, nfiles> file =
  {
    "Rootfiles/DeltaRPhi_cluster_upsilon_acts_full_no_distortion-new.root",
    "Rootfiles/DeltaRPhi_cluster_upsilon_acts_full_distorted-tony.root",
    "Rootfiles/DeltaRPhi_cluster_upsilon_acts_full_distorted_fullmap-tony.root"
  };

  // const TString pdfFile = "Figures/DeltaRPhi_compare_full_distortion-new.pdf";
  // const TString pdfFile = "Figures/DeltaRPhi_compare_truth_distortion-new.pdf";
  
  constexpr std::array<int, 3> color = { 1, 2, 4 };
  constexpr std::array<int, 3> symbol = { 20, 20, 20 };

  std::array<TString, nfiles> labels =
  {
    "no distortions",
    "beam-induced distortions + correction",
    "static distortions + correction"
  };
  
  const TString pdfFile = "Figures/DeltaRPhi_cluster_compare_full_distortions.pdf";

  PdfDocument pdfDocument( pdfFile );

  {
    const TString cvName( "cv" );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    cv->SetTopMargin( 0.07 );
    cv->SetLeftMargin( 0.16 );
    
    // auto h = new TH1F( "dummy", "", 100, 0, 85);
    auto h = new TH1F( "dummy", "", 100, 20, 80 );
    h->SetMinimum(0);
    h->SetMaximum(250);
    h->GetXaxis()->SetTitle( "r (cm)" );
    // h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-cluster) (#mum)" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.5 );
    h->GetYaxis()->SetMaxDigits(4);
    h->Draw();
    
    auto legend = new TLegend( 0.2, 0.18, 0.87, 0.35, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();
        
    for( int i = 0; i < nfiles; ++i )
    {
      std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
      auto tg = static_cast<TGraphErrors*>( input->Get( "residuals" ) );
      tg->SetMarkerStyle( symbol[i] );
      tg->SetMarkerColor( color[i] );
      tg->SetLineColor( color[i] );
      
      tg->Draw( "P" );
      
      legend->AddEntry( tg, labels[i], "AP" );      
    }
    
    pdfDocument.Add( cv );
  
  }
  
  {
    const TString cvName( "cv2" );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    cv->SetTopMargin( 0.07 );
    cv->SetLeftMargin( 0.16 );
    
    // auto h = new TH1F( "dummy", "", 100, 0, 85);
    auto h = new TH1F( "dummy", "", 100, 20, 80 );
    h->SetMinimum(-0.01);
    h->SetMaximum(0.01);
    h->GetXaxis()->SetTitle( "r (cm)" );
    // h->GetYaxis()->SetTitle( "#LTr.#Delta#phi#GT (track-cluster) (cm)" );
    h->GetYaxis()->SetTitle( "#LTr.#Delta#phi#GT (cluster-truth) (cm)" );
    h->GetYaxis()->SetTitleOffset( 1.5 );
    h->GetYaxis()->SetMaxDigits(4);
    h->Draw();
    
    auto legend = new TLegend( 0.2, 0.18, 0.87, 0.35, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();
        
    for( int i = 0; i < nfiles; ++i )
    {
      std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
      auto tg = static_cast<TGraphErrors*>( input->Get( "residuals mean" ) );
      tg->SetMarkerStyle( symbol[i] );
      tg->SetMarkerColor( color[i] );
      tg->SetLineColor( color[i] );
      
      tg->Draw( "P" );
      
      legend->AddEntry( tg, labels[i], "AP" );      
    }
    
    pdfDocument.Add( cv );
  
  }
}
