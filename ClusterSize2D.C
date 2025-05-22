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
TString ClusterSize2D( TString tag = TString() )
{

  set_style( false );

  // input files
  if( tag.IsNull() ) tag = "_2k_realistic_full_zhengyun" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/ClusterSize2D%s.pdf", tag.Data() );

  std::cout << "ClusterSize2D - inputFile: " << inputFile << std::endl;
  std::cout << "ClusterSize2D - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._z_size:_clusters._phi_size" );
  const TString var3d = Form( "%s:_clusters._layer", var.Data() );
  const TCut cluster_cut;
  const TCut momentum_cut;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // cluster size not available (yet) for TPC
    if( idet == 2 ) continue;

    const TString hname( Form( "ClusterSize2D_%i", idet ) );
    std::unique_ptr<TH3> h3d( new TH3F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 5, 0, 5, 5, 0, 5 ) );
    Utils::TreeToHisto( tree, hname, var3d, cluster_cut&&momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, (idet == 0 ) ? 400:800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      // auto h = h3d->ProjectionY( hname, ilayer+1, ilayer+1 );
      h3d->GetXaxis()->SetRange( ilayer+1, ilayer+1 );
      auto h = h3d->Project3D( "zy" );
      h->SetName( hname );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "csize_{#phi}" );
      h->GetYaxis()->SetTitle( "csize_{z}" );

      cv->cd( ilayer+1 );
      h->Draw( "colz" );
      gPad->SetRightMargin(0.2);
      gPad->Update();
      // gPad->SetLogy( true );
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
