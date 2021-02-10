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


  #if false
  constexpr int nFiles = 9;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi20k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nominal.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi5k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi2k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi1k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi500.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi200.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi100.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_nphi50.root"
  };

  const TString pdfFile = "Figures/DeltaRPhi_truth_realistic_truth_notpc_nphi.pdf";
  #else

  constexpr int nFiles = 9;
  const std::array<TString, nFiles> file =
  {
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi20k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nominal.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi5k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi2k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi1k.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi500.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi200.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi100.root",
    "Rootfiles/DeltaRPhi_truth_realistic_truth_notpc_single_nphi50.root"
  };

  const TString pdfFile = "Figures/DeltaRPhi_truth_realistic_truth_notpc_single_nphi.pdf";
  #endif

  constexpr std::array<int, nFiles> color = { kMagenta+1, kBlue+2, kBlue, kCyan+2, kGreen+2, kGreen+1, kYellow-4, kOrange-3, kRed };
  constexpr std::array<int, nFiles> symbol = { 20, 20, 20, 20, 20, 20, 20, 20, 20 };
  const std::array<TString, nFiles> label = 
  {
    "N_{r#phi}=20k, #sigma^{OT}_{r#phi}=40 #mum",
    "N_{r#phi}=10k, #sigma^{OT}_{r#phi}=110 #mum",
    "N_{r#phi}=5k, #sigma^{OT}_{r#phi}= 260 #mum",
    "N_{r#phi}=2k, #sigma^{OT}_{r#phi}=710 #mum",
    "N_{r#phi}=1k, #sigma^{OT}_{r#phi}=1.5 mm",
    "N_{r#phi}=500, #sigma^{OT}_{r#phi}=3 mm",
    "N_{r#phi}=200, #sigma^{OT}_{r#phi}=7.5 mm",
    "N_{r#phi}=100, #sigma^{OT}_{r#phi}=15 mm",
    "N_{r#phi}=50, #sigma^{OT}_{r#phi}=29 mm"
  };

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.16 );

  // auto h = new TH1F( "dummy", "", 100, 0, 90 );
  auto h = new TH1F( "dummy", "", 100, 20, 80 );
  h->SetMinimum(0);
  h->SetMaximum(12e3);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.5 );
  h->GetYaxis()->SetMaxDigits(4);
  h->Draw();

  auto legend = new TLegend( 0.16, 0.45, 0.60, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

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

    std::cout << "DrawDeltaRPhi - file: " << file[i] << " maximum: " << maximum << std::endl;

    maximumList.push_back( maximum );

  }

  // print maximum values
  Stream::PrintVector( "double", "residuals", maximumList, "%.0f" );

  cv->SaveAs( pdfFile );

}
