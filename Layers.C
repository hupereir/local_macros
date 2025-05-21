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
  // if( tag.IsNull() ) tag = "_genfit_truth";
  if( tag.IsNull() ) tag = "_acts_truth-new";
  // if( tag.IsNull() ) tag = "_acts_truth_realistic";
  // if( tag.IsNull() ) tag = "_acts_truth_realistic_notpc";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst*_1?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/Layers%s.pdf", tag.Data() );
  std::cout << "Layer - inputFile: " << inputFile << std::endl;
  std::cout << "Layer - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TCut trk_cut;
//   const TCut trk_cut( 
//     "_tracks._nclusters_tpc > 0 "
//     // "&& _tracks._nclusters_micromegas >= 1 "
//     "&& _tracks._is_primary "
//     "&& fabs(_tracks._eta)<1 "
//     "&& _tracks._px != 0 && _tracks._py !=0" );
//   // const TCut trk_cut( "_tracks._truth_pt >= 3 && _tracks._is_primary && _tracks._nclusters_tpc > 0" );

  {
    const TString var( "_tracks._clusters._layer" );
    const TString hname( "layer" );
    auto h( new TH1F( hname, "", nLayersTotal+1, 0, nLayersTotal+1 ) );
    Utils::TreeToHisto( tree, hname, var, trk_cut, false );
    
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( "layer" );
    
    std::cout << "Layers - entries: " << h->GetEntries() << std::endl;
    
    // normalization
    auto hnorm = Utils::TreeToHisto( tree, "norm", "_tracks._truth_pt", trk_cut, true );
    auto norm = hnorm->GetEntries();
    
    std::cout << "Layers - normalisation: " << norm << std::endl;
    h->Scale( 1./norm );
    h->SetMinimum(0);
    
    // create canvas
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    h->SetMarkerStyle( 20 );
    h->SetMarkerColor( 2);
    h->SetLineColor(2);
    h->Draw();
    
    Draw::HorizontalLine( cv, 1.0 )->Draw();

    Draw::PutText( 0.6, 0.6, Form( "N_{tracks} = %.0f", norm ) );
    
    pdfDocument.Add( cv );  
  }
  
  {
    
    auto h = new TH1F( "clusters_tpc", "", 50, 0, 50 );
    const TString var( "_tracks._nclusters_tpc" );
    Utils::TreeToHisto( tree, h->GetName(), var, trk_cut, false );
    h->GetXaxis()->SetTitle( "N_{clusters, tpc}" );
    
    auto cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    h->Draw();
    
    std::cout << "Layers - mean: " << h->GetMean() << std::endl;
    
    pdfDocument.Add( cv );  
  }
    
  return pdfFile;

}
