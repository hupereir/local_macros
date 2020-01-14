R__LOAD_LIBRARY(libRootUtilBase.so)
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>

void Layers( TString file = "rootfiles/G4sPHENIX_g4svtx_eval.root" )
{

  // layer histogram
  TH1* hLayer = new TH1F( "hLayer", "hLayer", 60, 0, 60 );
  auto chain = FileManager( file ).GetChain( "ntp_cluster" );
  Utils::TreeToHisto( chain, "hLayer", "layer", "", false );

  PdfDocument pdf( "Figures/Layers.pdf" );

  TCanvas* cv = new TCanvas( "cv", "cv", 700, 700 );
  hLayer->GetXaxis()->SetTitle( "Layer Id" );
  hLayer->Draw();

  pdf.Add( cv );

}
