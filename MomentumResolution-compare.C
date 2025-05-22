#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString MomentumResolution( TString tag = TString() )
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 2;
  const std::array<TString, 2> file =
  {
    "DST/dst_eval.root",
    "DST/CONDOR_full_high_occupancy_new/dst_eval_full_high_occupancy_new*.root"
  };

  const TString pdfFile = "Figures/MomentumResolution_compare.pdf";

  constexpr std::array<int, 2> color = { 1, 2 };
  constexpr std::array<int, 2> symbol = { 20, 20 };
  const std::array<TString, 2> label =
  {
    "old errors",
    "new errors"
  };

  // variable names
  const TString var( "_tracks._pt/_tracks._truth_pt" );
  const TCut cut( "abs(_tracks._truth_pt)>0.5&&_tracks._nclusters_mvtx==3&&_tracks._nclusters_intt>=2" );

  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  cv->SetLeftMargin( 0.16 );

  auto legend = new TLegend( 0.5, 0.84, 0.97, 0.94, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  for( int i = 0; i < nFiles; ++i )
  {
    std::cout << "MomentumResolution - processing " << file[i] << std::endl;
    
    FileManager fileManager( file[i] );
    auto tree = fileManager.GetChain( "T" );
    
    gROOT->cd();
    const auto hname = Form( "h%i", i );
    auto h = new TH1F( hname, "", 100, 0.7, 1.3 );
    Utils::TreeToHisto( tree, hname, var, cut, false );
    h->SetTitle( "" );
    h->SetMaximum( 4e3 );
    h->GetXaxis()->SetTitle( "#it{p}_{T,track}/#it{p}_{T,truth}" );

    h->SetLineColor( color[i] );
    if( i ) h->Draw( "same" );
    else h->Draw();

    std::cout << "MomentumResolution - RMS: " << h->GetRMS() << " +/- " << h->GetRMSError() << std::endl;
    legend->AddEntry( h, Form( "%s, RMS=%.4f#pm%.4f", label[i].Data(), h->GetRMS(), h->GetRMSError() ), "PL" );
    legend->Draw();

    gPad->SetLogy( true );

  }

  cv->SaveAs( pdfFile );
  return pdfFile;
}
