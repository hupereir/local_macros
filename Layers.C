#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString Layers( TString tag = TString() )
{

  // if( tag.IsNull() ) tag = "_full";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/Layers%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/Layers%s.root", tag.Data() );

  std::cout << "Layer - inputFile: " << inputFile << std::endl;
  std::cout << "Layer - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._layer" );
  const TCut trk_cut( "_tracks._truth_pt >= 3 && _tracks._is_primary" );

  const TString hname( "layer" );
  std::unique_ptr<TH1> h( new TH1F( hname, "", nLayersTotal, 0, nLayersTotal ) );
  Utils::TreeToHisto( tree, hname, var, trk_cut, false );

  h->SetTitle( "" );
  h->GetXaxis()->SetTitle( "layer" );

  // normalization
  auto hnorm = Utils::TreeToHisto( tree, "norm", "_tracks._truth_pt", trk_cut, true );
  auto norm = hnorm->GetEntries();
  h->Scale( 1./norm );
  h->SetMinimum(0);
  
  // create canvas
  const TString cvName = "cv";
  std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
  h->Draw();

  Draw::HorizontalLine( cv.get(), 1.0 )->Draw();
  
  pdfDocument.Add( cv.get() );

  // save everything in rootfiles
  {
    std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
    output->cd();
    h->Write();
    output->Close();
  }

  return pdfFile;

}
