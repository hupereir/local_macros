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
void DrawDeltaZ()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 5;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaZ_truth_flat_full_notpc_single_nz11k_highpt.root",
    "Rootfiles/DeltaZ_truth_flat_full_notpc_single_nominal_highpt.root",
    "Rootfiles/DeltaZ_truth_flat_full_notpc_single_nz3k_highpt.root",
    "Rootfiles/DeltaZ_truth_flat_full_notpc_single_nz1k_highpt.root",
    "Rootfiles/DeltaZ_truth_flat_full_notpc_single_nz500_highpt.root"
  };

  const TString pdfFile = "Figures/DeltaZ_truth_notpc_single_nz_highpt-flat.pdf";

//
//   constexpr int nFiles = 5;
//   const std::array<TString, nFiles> file =
//   {
//     "Rootfiles/DeltaZ_truth_flat_full_notpc_nz11k_highpt.root",
//     "Rootfiles/DeltaZ_truth_flat_full_notpc_nominal_highpt.root",
//     "Rootfiles/DeltaZ_truth_flat_full_notpc_nz3k_highpt.root",
//     "Rootfiles/DeltaZ_truth_flat_full_notpc_nz1k_highpt.root",
//     "Rootfiles/DeltaZ_truth_flat_full_notpc_nz500_highpt.root"
//   };
//
//   const TString pdfFile = "Figures/DeltaZ_truth_notpc_nz_highpt-flat.pdf";

  PdfDocument pdfDocument( pdfFile );

  constexpr std::array<int, 5> color = { kBlue, kCyan+2, kGreen+2, kOrange+1, kRed };
  constexpr std::array<int, 5> symbol = { 20, 20, 20, 20, 20 };
  const std::array<TString, 5> label =
  {
    "n_{z}=11k, #sigma_{z}=40 #mum",
    "n_{z}=5k, #sigma_{z}=100 #mum",
    "n_{z}=3k, #sigma_{z}= 220 #mum",
    "n_{z}=1k, #sigma_{z}=560 #mum",
    "n_{z}=500, #sigma_{z}=1.2 mm"
  };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.15 );

  // auto h = new TH1F( "dummy", "", 100, 0, 90 );
  auto h = new TH1F( "dummy", "", 100, 20, 80 );
  h->SetMinimum(0);
  // h->SetMaximum(12.5e2);
  h->SetMaximum(5.5e2);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.4 );
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
