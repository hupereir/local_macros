#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

TString DrawMomentumResolution_2d()
{
  gStyle->SetOptStat(0);
  set_style( false );

  static constexpr int nfiles = 2;
  std::array<TString, nfiles> filenames =
  {
    "Rootfiles/MomentumResolution_2d_flat_acts_truth_nodistortion.root",
    "Rootfiles/MomentumResolution_2d_flat_acts_truth_notpc_nodistortion_smeared_mm.root"
  };

  std::array<TString, nfiles> labels =
  {
    "MVTX+INTT+TPC",
    "MVTX+INTT+TPOT"
  };

  std::array<int, nfiles> colors = {1, 2};
  const TString pdfFile = "Figures/MomentumResolution_2d_flat_truth_nodistortion-compare_smeared.pdf";

  PdfDocument pdfDocument( pdfFile );

  {
    // resolution
    TCanvas* cv = new TCanvas( "cv", "cv", 800, 600 );
    cv->SetLeftMargin( 0.17 );

    auto h = new TH1F( "dummy", "", 100, 0, 50 );
    h->SetMinimum(0);
    h->SetMaximum(0.1);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "#sigma( #it{p}_{T,track}/#it{p}_{T,truth} )" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    auto legend = new TLegend( 0.2, 0.7, 0.87, 0.85, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    for( int i = 0; i < nfiles; ++i )
    {
      auto file = TFile::Open( filenames[i] );
      auto tg = static_cast<TGraphErrors*>(  file->Get("resolution") );
      tg->SetMarkerStyle( 20 );
      tg->SetMarkerColor( colors[i] );
      tg->SetLineColor( colors[i] );
      tg->Draw("P" );
      legend->AddEntry( tg, labels[i], "PL" );
    }

    pdfDocument.Add( cv );
  }

  {
    // momentum scale
    TCanvas* cv = new TCanvas( "cv2", "cv2", 800, 600 );
    cv->SetLeftMargin( 0.17 );

    auto h = new TH1F( "dummy", "", 100, 0, 20 );
    h->SetMinimum(0.8);
    h->SetMaximum(1.1);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "#LT#it{p}_{T,track}/#it{p}_{T,truth}#GT" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    auto legend = new TLegend( 0.2, 0.18, 0.87, 0.35, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    for( int i = 0; i < nfiles; ++i )
    {
      auto file = TFile::Open( filenames[i] );
      auto tg = static_cast<TGraphErrors*>(  file->Get("momentum scale") );
      tg->SetMarkerStyle( 20 );
      tg->SetMarkerColor( colors[i] );
      tg->SetLineColor( colors[i] );
      tg->Draw("P" );
      legend->AddEntry( tg, labels[i], "PL" );
    }

    Draw::HorizontalLine( cv, 1.0 )->Draw();

    pdfDocument.Add( cv );
  }

  return pdfFile;
}
