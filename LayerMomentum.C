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
TString LayerMomentum TString tag = TString() )
{

  const auto maxMomentum = 5;

  if( tag.IsNull() ) tag = "_truth_notpc_noouter";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/Momentum_truth%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/Momentum_truth%s.root", tag.Data() );

  std::cout << "LayerMomentum - inputFile: " << inputFile << std::endl;
  std::cout << "LayerMomentum - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._truth_pt" );
  const TString var2d = Form( "%s:_tracks._clusters[]._layer", var.Data() );
  const TCut momentum_cut = Form( "_tracks._truth_pt < %d", maxMomentum );

  const TString hname( "momentum" );
  std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayersTotal, 0, nLayersTotal, 100, 0, maxMomentum ) );
  Utils::TreeToHisto( tree, hname, var2d, momentum_cut, false );

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;

      const auto hname = Form( "h_%i", layerIndex );
      std::unique_ptr<TH1> h( h2d->ProjectionY( hname, layerIndex+1, layerIndex+1 ) );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );

      cv->cd( ilayer+1 );
      h->Draw();

      Draw::PutText( 0.2, 0.8, Form( "layer = %i, entries: %.0f", layerIndex, h->GetEntries() ) );
      Draw::PutText( 0.2, 0.7, Form( "<#it{p}_{T}> = %.3f GeV/#it{c}", h->GetMean() ) );

      gPad->Update();

      // save in array
      h_array[layerIndex] = std::move(h);

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  // save everything in rootfiles
  {
    std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
    output->cd();
    for( auto&& h:h_array) { h->Write(); }
    output->Close();
  }

  return pdfFile;

}
