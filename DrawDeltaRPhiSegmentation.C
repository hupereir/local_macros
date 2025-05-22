#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
void DrawDeltaRPhiSegmentation()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   constexpr int nFiles = 9;
//   const std::array<TString, nFiles> file =
//   {
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi20k.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nominal.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi5k.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi2k.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi1k.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi500.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi200.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi100.root",
//     "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi50.root"
//   };
//
//   const TString pdfFile = "Figures/DeltaRPhiSegmentation_truth_realistic_truth_notpc_nphi.pdf";

  constexpr int nFiles = 9;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi20k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nominal.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi5k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi2k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi1k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi500.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi200.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi100.root",
    "Rootfiles/DeltaRPhi_truth_realistic_full_notpc_nphi50.root"
  };

  const TString pdfFile = "Figures/DeltaRPhiSegmentation_truth_realistic_full_notpc_nphi.pdf";

  constexpr std::array<double, nFiles> segmentation = { 2e4, 1e4, 5e3, 2e3, 1e3, 500, 200, 100, 50 };

  PdfDocument pdfDocument( pdfFile );

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.16 );

  // auto h = new TH1F( "dummy", "", 100, 0, 90 );
  auto h = new TH1F( "dummy", "", 100, 10, 3e5 );
  h->SetMinimum(0);
  h->SetMaximum(1e4);
  h->GetXaxis()->SetTitle( "N_{r#phi}" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.5 );
  h->GetYaxis()->SetMaxDigits(4);
  h->Draw();
  gPad->SetLogx(true);

  auto tgloc = new TGraphErrors();
  for( int i = 0; i < nFiles; ++i )
  {
    std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
    auto tg = static_cast<TGraphErrors*>( input->Get( "residuals" ) );
    tgloc->SetPoint( i, segmentation[i], Utils::GetMaximum( tg ) );
    tgloc->SetPointError( i, 0, Utils::GetMaximumError( tg ) );
  }

  tgloc->SetMarkerStyle( 20 );
  tgloc->SetMarkerColor( 1 );
  tgloc->Draw( "P" );
  pdfDocument.Add( cv );

}
