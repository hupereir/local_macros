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
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

//____________________________________________________________________________
TString ClusterSize( TString tag = TString() )
{

  set_style( false );

//   // input files
//   std::array<int, nDetectors> max_csize = { 20, 10, 10, 10, 10, 10, 10 };
//   if( tag.IsNull() ) tag = "_1k_realistic_micromegas";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  // input files
//   std::array<int, nDetectors> max_csize = { 20, 10, 10, 10, 10, 30, 30 };
//   if( tag.IsNull() ) tag = "_Hijing_micromegas_merged";
//   const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_merged/dst_eval_sHijing_0-12fm_merged_0*root";

//   // input files
//   std::array<int, nDetectors> max_csize = { 20, 10, 10, 10, 10, 30, 30 };
//   if( tag.IsNull() ) tag = "_Hijing_micromegas_single";
//   const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval/dst_eval_sHijing_0-12fm_0*root";
  
  std::array<int, nDetectors> max_csize = { 20, 10, 10, 10, 10, 30, 30 };
  if( tag.IsNull() ) tag = "_genfit_truth_realistic-newgeom";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/ClusterSize%s.pdf", tag.Data() );
  const TString rootFile = Form( "Rootfiles/ClusterSize%s.root", tag.Data() );
  
  std::cout << "ClusterSize - inputFile: " << inputFile << std::endl;
  std::cout << "ClusterSize - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._size" );
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );

  const TString var_tpc( "_clusters._phi_size" );
  const TString var2d_tpc = Form( "%s:_clusters._layer", var_tpc.Data() );

  const TCut cluster_cut;
  const TCut momentum_cut;

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> histograms;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    const TString hname( Form( "clusterSize_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], max_csize[idet], 0, max_csize[idet] ) );
    Utils::TreeToHisto( tree, hname, is_tpc(idet) ? var2d_tpc: var2d, cluster_cut&&momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv;
    if( idet == 0 ) cv.reset( new TCanvas( cvName, cvName, 900, 300 ) );
    else if( idet >= 5 ) cv.reset( new TCanvas( cvName, cvName, 400, 400 ) );
    else  cv.reset( new TCanvas( cvName, cvName, 900, 900 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      std::unique_ptr<TH1> h( h2d->ProjectionY( hname, ilayer+1, ilayer+1 ) );
      h->SetTitle( "" );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( is_tpc(idet) ? Form( "csize_{#phi} (layer %i)", layerIndex ) : Form( "csize (layer %i)", layerIndex ) );
      h->SetMaximum( h->GetMaximum()*1.2 );

      cv->cd( ilayer+1 );
      h->Draw();
      
      auto mean = h->GetMean();
      auto meanError = h->GetMeanError();
      Draw::VerticalLine( gPad, mean )->Draw();
      Draw::PutText( 0.5, 0.7, Form( "Mean = %.1f", mean ) );
      gPad->SetLogy(true);
      gPad->Update();    
    
      // save in array
      histograms[layerIndex] = std::move(h);
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }
  
  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:histograms) { if(h) h->Write(); }
  output->Close();

  return pdfFile;
}
