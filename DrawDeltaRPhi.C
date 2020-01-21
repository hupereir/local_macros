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
void DrawDeltaRPhi()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 5;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaRPhi_truth_flat_full_notpc_single_nphi20k_highpt.root",
    "Rootfiles/DeltaRPhi_truth_flat_full_notpc_single_nominal_highpt.root",
    "Rootfiles/DeltaRPhi_truth_flat_full_notpc_single_nphi5k_highpt.root",
    "Rootfiles/DeltaRPhi_truth_flat_full_notpc_single_nphi2k_highpt.root",
    "Rootfiles/DeltaRPhi_truth_flat_full_notpc_single_nphi1k_highpt.root"
  };

  const TString pdfFile = "Figures/DeltaRPhi_truth_notpc_single_nphi_highpt-flat.pdf";
  PdfDocument pdfDocument( pdfFile );

//   constexpr int nFiles = 5;
//   const std::array<TString, nFiles> file =
//   {
//     "Rootfiles/DeltaRPhi_truth_flat_full_notpc_nphi20k_highpt.root",
//     "Rootfiles/DeltaRPhi_truth_flat_full_notpc_nominal_highpt.root",
//     "Rootfiles/DeltaRPhi_truth_flat_full_notpc_nphi5k_highpt.root",
//     "Rootfiles/DeltaRPhi_truth_flat_full_notpc_nphi2k_highpt.root",
//     "Rootfiles/DeltaRPhi_truth_flat_full_notpc_nphi1k_highpt.root"
//   };
//
//   const TString pdfFile = "Figures/DeltaRPhi_truth_notpc_nphi_highpt-flat.pdf";
//   PdfDocument pdfDocument( pdfFile );

  constexpr std::array<int, 5> color = { kBlue, kCyan+2, kGreen+2, kOrange+1, kRed };
  constexpr std::array<int, 5> symbol = { 20, 20, 20, 20, 20 };
  const std::array<TString, 5> label = {
    "n_{#phi}=20k, #sigma_{r#phi}=40 #mum",
    "n_{#phi}=10k, #sigma_{r#phi}=110 #mum",
    "n_{#phi}=5k, #sigma_{r#phi}= 260 #mum",
    "n_{#phi}=2k, #sigma_{r#phi}=710 #mum",
    "n_{#phi}=1k, #sigma_{r#phi}=1.5 mm"
  };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.16 );

  // auto h = new TH1F( "dummy", "", 100, 0, 90 );
  auto h = new TH1F( "dummy", "", 100, 20, 80 );
  h->SetMinimum(0);
  h->SetMaximum(1.5e3);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.5 );
  h->Draw();

  auto legend = new TLegend( 0.16, 0.6, 0.60, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  for( int i = 0; i < nFiles; ++i )
  {
    std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
    auto tg = static_cast<TGraphErrors*>( input->Get( "residuals" ) );
    tg->SetMarkerStyle( symbol[i] );
    tg->SetMarkerColor( color[i] );
    tg->SetLineColor( color[i] );

    tg->Draw( "P" );

    legend->AddEntry( tg, label[i], "AP" );

  }

  pdfDocument.Add( cv );
}
