#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

void TrackingEfficiency_QA()
{
  set_style( false );
  gStyle->SetOptStat(0);


  const TString flag = "_low_multiplicity-new3";
  const TString inputFile = Form( "Rootfiles/QA%s/TrackingQA_*.root", flag.Data());
  const TString pdfFile = Form( "Figures/TrackingEfficiency_QA%s.pdf", flag.Data());

  std::cout << "TrackingEfficiency_QA - inputFile: " << inputFile << std::endl;
  std::cout << "TrackingEfficiency_QA - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  FileManager fileManager( inputFile );

  const TString href_name = "h_QAG4SimulationTracking_nGen_pTGen";
  const TString hfound_name = "h_QAG4SimulationTracking_nMVTX_nReco_pTGen";

  const auto href = fileManager.GetHistogram(href_name);
  const auto hfound = static_cast<TH2*>(fileManager.GetHistogram(hfound_name));
  auto heff = static_cast<TH1*>( href->Clone("href" ));
  heff->Reset();

  // require at least 2 MVTX cluster
  const int min_mvtx_clusters = 2;
  const int bin_start = hfound->GetYaxis()->FindBin( min_mvtx_clusters );
  const auto hfound_1d = hfound->ProjectionX(Form( "%s_proj", hfound->GetTitle()), bin_start);

  Utils::DivideHistograms(hfound_1d, href, heff);
  heff->SetMarkerStyle(20);
  heff->GetXaxis()->SetTitle( "p_{T} (GeV/c)" );
  heff->GetYaxis()->SetTitle( "efficiency" );
  heff->SetMaximum(1.1);
  heff->SetMinimum(0);

  if( false )
  {
    const double ptmin = 0;
    const double ptmax = 10;
    heff->GetXaxis()->SetRangeUser( ptmin, ptmax );
  }

  auto cv = new TCanvas( "cv", "cv", 900, 900 );
  gPad->SetTopMargin( 0.07 );
  gPad->SetLeftMargin( 0.15);
  gPad->SetRightMargin( 0.05);

  heff->Draw();

  gPad->SetLogx(true);
  gPad->SetGridy(true);
  gPad->Update();
  Draw::HorizontalLine( gPad, 1 )->Draw();

  pdfDocument.Add(cv);
}
