#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString HitStatistics()
{

  set_style( false );

  // const TString tag = "_upsilon_acts_truth_distorted_fullmap";
  // const TString tag = "_upsilon_acts_full_no_distortion";
  // const TString tag = "_upsilon_acts_truth_distorted_test";
  // const TString tag = "_upsilon_acts_full_no_distortion_nomm";
  const TString tag = "_upsilon_acts_full_no_distortion_thinmm";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/HitStatistics%s.pdf", tag.Data() );

  std::cout << "HitStatistics - inputFile: " << inputFile << std::endl;
  std::cout << "HitStatistics - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) { 
    std::cout << "HitStatistics - invalid tree" << std::endl;
    return pdfFile;
  }

  // variable names
  {
    const TString var( "_events._nclusters_tpc" );
    const TCut cut = Form( "%s< 1000", var.Data() );
    
    auto h = Utils::TreeToHisto( tree, "nclusters_tpc", var, cut, true );
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    h->Draw();
    pdfDocument.Add( cv );
    
    std::cout << "HitStatistics - " << var << " mean: " << h->GetMean() << " +/- " << h->GetMeanError() << std::endl;
  }
  

  // variable names
  {
    const TString var( "_events._ng4hits_tpc" );
    const TCut cut = Form( "%s< 2000", var.Data() );
    
    auto h = Utils::TreeToHisto( tree, "nclusters_tpc", var, cut, true );
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    h->Draw();
    pdfDocument.Add( cv );
    
    std::cout << "HitStatistics - " << var << " mean: " << h->GetMean() << " +/- " << h->GetMeanError() << std::endl;
  }

  return pdfFile;
}
