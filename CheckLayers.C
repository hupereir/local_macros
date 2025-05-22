#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

TTree* tree = nullptr;

//_________________________________________
TCanvas* Draw( TString var1, TString var2, int max_layers )
{   
  
  TH1* h1 = new TH1F( "h1", "h1", max_layers, 0, max_layers );
  Utils::TreeToHisto( tree, h1->GetName(), var1, TCut(), false );
  
  TH1* h2 = new TH1F( "h2", "h2", max_layers, 0, max_layers );
  Utils::TreeToHisto( tree, h2->GetName(), var2, TCut(), false );
  
  auto cv = new TCanvas( "cv", "cv", 1200, 600 );
  cv->Divide( 2, 1 );
  
  cv->cd(1);
  h1->SetLineColor( 1 );
  h1->Draw();
  
  h2->SetLineColor(2);
  h2->Draw("same");
  
  cv->cd(2);
  const auto vardiff = Form( "(%s)-(%s)", var1.Data(), var2.Data() );
  TH1* h3 = new TH1F( "h3", "h3", 4, -2, 2 );
  Utils::TreeToHisto( tree, h3->GetName(), vardiff, TCut(), false );
  h3->Draw("same");
  
  return cv;
}

//_________________________________________
TString CheckLayers()
{
  
  const TString tag = "_flat_truth_acts";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/CheckLayers%s.pdf", tag.Data() );

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  tree = fileManager.GetChain( "T" );
  if( !tree ) { 
    std::cout << "CheckLayers - invalid tree" << std::endl;
    return pdfFile;
  }
  
  {
    const TString var1 = "_tracks._nclusters_mvtx";
    const TString var2 = "TrackingEvaluator_hp::get_nclusters_mvtx( _tracks._mask )";
    pdfDocument.Add( Draw( var1, var2, 5 ) );
  }
    
  {
    const TString var1 = "_tracks._nclusters_intt";
    const TString var2 = "TrackingEvaluator_hp::get_nclusters_intt( _tracks._mask )";
    pdfDocument.Add( Draw( var1, var2, 5 ) );
  }
    
  {
    const TString var1 = "_tracks._nclusters_tpc";
    const TString var2 = "TrackingEvaluator_hp::get_nclusters_tpc( _tracks._mask )";
    pdfDocument.Add( Draw( var1, var2, 50 ) );
  }

  {
    const TString var1 = "_tracks._nclusters_micromegas";
    const TString var2 = "TrackingEvaluator_hp::get_nclusters_micromegas( _tracks._mask )";
    pdfDocument.Add( Draw( var1, var2, 5 ) );
  }

  return pdfFile;

}
