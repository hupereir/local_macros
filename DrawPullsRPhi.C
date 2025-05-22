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
void DrawPullsRPhi()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 2;
//   const std::array<TString, nFiles> file =
//   {
//     "Rootfiles/PullsRPhi_cluster_2k_realistic_full_nominal_old.root",
//     "Rootfiles/PullsRPhi_cluster_2k_realistic_full_nominal_new.root"
//   };
//   const TString pdfFile = "Figures/PullsRPhi_compare_new.pdf";
//
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/PullsRPhi_truth_2k_realistic_full_nominal_old.root",
    "Rootfiles/PullsRPhi_truth_2k_realistic_full_nominal_new.root"
  };
  const TString pdfFile = "Figures/PullsRPhi_truth_compare_new.pdf";

  constexpr std::array<int, nFiles> color = { 1, 2 };
  constexpr std::array<int, nFiles> symbol = { 20, 20 };
  const std::array<TString, nFiles> label = {
    "without fix",
    "with fix"
  };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 600 );
  cv->SetLeftMargin( 0.16 );

  auto h = new TH1F( "dummy", "", 100, 0, 90 );
  h->SetMinimum(0);
  h->SetMaximum(2);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "pulls_{r.#Delta#phi} (cluster-truth)" );
  h->GetYaxis()->SetTitleOffset( 1.5 );
  h->GetYaxis()->SetMaxDigits(4);
  h->Draw();

  auto legend = new TLegend( 0.16, 0.7, 0.60, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  std::vector<double> maximumList;

  for( int i = 0; i < nFiles; ++i )
  {
    std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
    auto tg = static_cast<TGraphErrors*>( input->Get( "pulls" ) );
    tg->SetMarkerStyle( symbol[i] );
    tg->SetMarkerColor( color[i] );
    tg->SetLineColor( color[i] );

    tg->Draw( "P" );

    legend->AddEntry( tg, label[i], "AP" );

    // get the maximum in the tpc
    double maximum = 0;
    for( int ip = 0; ip < tg->GetN(); ++ip )
    {
      double r, sigma;
      tg->GetPoint( ip, r, sigma );
      if( r < rmin_tpc || r > rmax_tpc ) continue;
      if( sigma > maximum ) maximum = sigma;
    }

    std::cout << "DrawDeltaRPhi - file: " << file[i] << " maximum: " << maximum << std::endl;

    maximumList.push_back( maximum );

  }

  cv->SaveAs( pdfFile );

}
