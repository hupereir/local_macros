#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
TString TrackEfficiency_intt()
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  
  static constexpr int nFiles = 2;
  const std::array<TString, nFiles> file =
  {
    "DST/CONDOR_full_high_occupancy_old/dst_eval_full_high_occupancy_old*.root",
    "DST/CONDOR_full_high_occupancy_new/dst_eval_full_high_occupancy_new*.root"
  };

  const TString pdfFile = "Figures/TrackingEfficiency_intt_compare.pdf";

  static constexpr std::array<int, nFiles> color = { 1, 2 };
  static constexpr std::array<int, nFiles> symbol = { 20, 20 };
  const std::array<TString, nFiles> label =
  {
    "old errors",
    "new errors"
  };
 
  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );

  const TString var = "_tracks._pt";
  const TCut refCut = "_tracks._embed == 2";
  const TCut foundCut = "_tracks._nclusters_intt>=2";
  
  // loop over files
  for( int i = 0; i < nFiles; ++i )
  {
    std::cout << "TrackHits - processing " << file[i] << std::endl;
    
    FileManager fileManager( file[i] );
    auto tree = fileManager.GetChain( "T" );
    gROOT->cd();

    // reference histogram
    auto hname = Form( "href%i", i );
    auto href = new TH1F( hname, "", 40, &Utils::LogAxis( 40, 0.1, 50 )[0] );
    Utils::TreeToHisto( tree, hname, var, refCut, false );
    
    hname = Form( "hfound%i", i );
    auto hfound = new TH1F( hname, "", 40, &Utils::LogAxis( 40, 0.1, 50 )[0] );
    Utils::TreeToHisto( tree, hname, var, refCut&&foundCut, false );
    
    hname = Form( "heff%i", i );
    auto heff = new TH1F( hname, "", 40, &Utils::LogAxis( 40, 0.1, 50 )[0] );
    
    Utils::DivideHistograms( hfound, href, heff );

    
    heff->SetMarkerStyle( symbol[i] );
    heff->SetMarkerColor( color[i] );
    heff->SetLineColor( color[i] );
    
    if( i ) heff->Draw( "same" );
    else 
    {
      heff->Draw();
      gPad->SetLogx(true);
    }
    
  }
  
  cv->SaveAs( pdfFile );
  return pdfFile;
  
}
