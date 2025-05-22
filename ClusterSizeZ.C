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
TString ClusterSizeZ()
{

  set_style( false );

  // input files
  const TString tag = "_acts_full_notpc_nodistortion-new3";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco*_?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/ClusterSizeZ%s.pdf", tag.Data() );

  std::cout << "ClusterSizeZ - inputFile: " << inputFile << std::endl;
  std::cout << "ClusterSizeZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._z_size" );
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );
  const TCut cluster_cut;
  const TCut momentum_cut;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    const TString hname( Form( "clusterSizeZ_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 10, 0, 10 ) );
    Utils::TreeToHisto( tree, hname, var2d, cluster_cut&&momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      TH1* h = h2d->ProjectionY( hname, ilayer+1, ilayer+1 );
      h->SetTitle( hname );
      if( !h->GetEntries() ) continue;
       
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "csize_{z}" );
      h->SetMaximum( h->GetMaximum()*1.2 );

      cv->cd( ilayer+1 );
      h->Draw();
      gPad->Update();
      // gPad->SetLogy( true );
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
