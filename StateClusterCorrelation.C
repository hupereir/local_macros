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

#include <trackbase/TrkrClusterContainer.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


//____________________________________________________________________________
TString StateClusterCorrelation()
{

  set_style( false );

  // gStyle->SetOptStats(0);

  const TString tag = "_TrackReconstruction_genfit";
  const TString inputFile = Form( "DST/CONDOR%s/TRACKS-00053877-000?-full.root", tag.Data() );

  const TString pdfFile = Form( "Figures/StateClusterCorrelation%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/StateClusterCorrelation%s.root", tag.Data() );

  std::cout << "StateClusterCorrelation - inputFile: " << inputFile << std::endl;
  std::cout << "StateClusterCorrelation - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  Utils::max_entries = 100000;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const TString var( "_tracks._nclusters_tpc:_tracks._nstates_tpc" );
  const TCut cut = ("_tracks._nstates_tpc>0" );

  auto h = new TH2F( "h0", "", 48, 0, 48, 48, 0, 48 );
  Utils::TreeToHisto( tree, h->GetName(), var, cut, false );

  h->SetStats(0);

  auto cv( new TCanvas( "cv", "cv", 900, 900 ) );
  cv->SetRightMargin(0.2);
  h->Draw("colz");

  pdfDocument.Add(cv);
  return pdfFile;
}
