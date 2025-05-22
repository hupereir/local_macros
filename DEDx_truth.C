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


TString dedx_truth()
{

  set_style( false );
  gStyle->SetOptStat(0);

  // const TString tag = "_single_electron";
  // const TString tag = "_single_piminus";
  const TString tag = "_single_proton";

  const TString inputFile = Form( "DST/CONDOR%s/DST_RECO_lpt-new/dst_reco*.root", tag.Data() );
  const TString pdfFile = Form( "Figures/dedx_truth%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  std::cout << "dedx - entries: " << tree->GetEntries() << std::endl;
  Utils::max_entries = 10000;

  // const TString var( "_tracks._dedx:_tracks._p" );
  const TString var( "_tracks._truth_dedx:_tracks._p" );
  const TCut cut = "_tracks._p < 5";
  // const TCut pid_cut = "_tracks._pid == 11";
  const TCut pid_cut = "_tracks._pid == 2212";

  auto h = new TH2F( "dedx", "", 100, 0, 5, 100, 0, 10 );
  h->GetXaxis()->SetTitle( "q.p [GeV]" );
  h->GetYaxis()->SetTitle( "dE/dx [ADC/cm]" );
  h->GetYaxis()->SetTitleOffset(1.6);

  Utils::TreeToHisto( tree, h->GetName(), var, cut && pid_cut, false );

  auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.15);

  h->Draw( "colz" );
  pdfDocument.Add(cv);

  return pdfFile;

}
