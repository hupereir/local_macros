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
void DrawDeltaZ()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  #if false
  constexpr int nFiles = 9;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz11k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nominal.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz3k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz1k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz500.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz300.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz100.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz50.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_nz10.root"
  };

  const TString pdfFile = "Figures/DeltaZ_truth_realistic_truth_notpc_nz.pdf";

  #else

  constexpr int nFiles = 9;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz11k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nominal.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz3k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz1k.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz500.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz300.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz100.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz50.root",
    "Rootfiles/DeltaZ_truth_realistic_truth_notpc_single_nz10.root"
  };

  const TString pdfFile = "Figures/DeltaZ_truth_realistic_truth_notpc_single_nz.pdf";

  #endif

  constexpr std::array<int, nFiles> color = { kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, kBlue, kCyan+2, kGreen+2, kOrange+1 };
  constexpr std::array<int, nFiles> symbol = { 20, 20, 20, 20, 20, 21, 21, 21, 21 };
  const std::array<TString, nFiles> label =
  {
    "N_{z}=11k, #sigma^{OT}_{z}=40 #mum",
    "N_{z}=5k, #sigma^{OT}_{z}=100 #mum",
    "N_{z}=3k, #sigma^{OT}_{z}= 220 #mum",
    "N_{z}=1k, #sigma^{OT}_{z}=560 #mum",
    "N_{z}=500, #sigma^{OT}_{z}=1.2 mm",
    "N_{z}=300, #sigma^{OT}_{z}=2.3 mm",
    "N_{z}=100, #sigma^{OT}_{z}=6 mm",
    "N_{z}=50, #sigma^{OT}_{z}=12 mm",
    "N_{z}=10, #sigma^{OT}_{z}=60 mm"
  };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.15 );

  // auto h = new TH1F( "dummy", "", 100, 0, 90 );
  auto h = new TH1F( "dummy", "", 100, 20, 80 );
  h->SetMinimum(0);
  h->SetMaximum(5000);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.4 );
  h->GetYaxis()->SetMaxDigits(3);
  h->Draw();

  auto legend = new TLegend( 0.16, 0.48, 0.60, 0.93, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  std::vector<double> maximumList;

  for( int i = 0; i < nFiles; ++i )
  {
    std::unique_ptr<TFile> input( TFile::Open( file[i] ) );
    auto tg = static_cast<TGraphErrors*>( input->Get( "residuals" ) );
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

    std::cout << "DrawDeltaZ - file: " << file[i] << " maximum: " << maximum << std::endl;

    maximumList.push_back( maximum );

  }

  legend->Draw();

  // print maximum values
  Stream::PrintVector( "double", "residuals", maximumList, "%.0f" );

  cv->SaveAs( pdfFile );

}
