#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
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
TString Radius( TString tag = TString() )
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // pdf output
  // if( tag.IsNull() ) tag = "_realistic_truth_notpc_noouter" ;
  if( tag.IsNull() ) tag = "_flat_full_notpc_nominal_highpt" ;
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/Radius%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // configuration
  const Bool_t doFit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_r" );
  const TString var2d = Form( "%s:_layer", var.Data() );

  const float maxRadius = 90;
  const TString hName( "radius" );
  std::unique_ptr<TH2> h2d( new TH2F( hName, "", nLayersTotal, 0, nLayersTotal, maxRadius*100, 0, maxRadius ) );
  Utils::TreeToHisto( tree, hName, var2d, TCut(), kFALSE );

  std::vector<float> layerRadius;

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();

  // loop over detectors
  for( Int_t idet = 0; idet < nDetectors; ++idet )
  {

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    Draw::DivideCanvas( cv, nLayers[idet], kFALSE );

    // loop over layers
    for( Int_t ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      Int_t layerIndex = firstLayer[idet] + ilayer;
      const auto hName = Form( "h_%i", layerIndex );
      TH1* h = h2d->ProjectionY( hName, layerIndex+1, layerIndex+1 );
      h->SetTitle( hName );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "r (cm)" );

      cv->cd( ilayer+1 );
      h->Draw();

      auto radius = h->GetMean();
      auto radius_error = h->GetMeanError();
      tg->SetPoint( layerIndex, layerIndex, radius );
      tg->SetPointError( layerIndex, 0, radius_error );

      layerRadius.push_back( radius );

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv );

  }

  // TGraph
  auto cv = new TCanvas( "cvtg", "cvtg", 800, 800 );
  cv->SetLeftMargin( 0.15 );

  auto h = new TH1F( "dummy", "", nLayersTotal, 0, nLayersTotal );
  h->SetMinimum(0);
  h->SetMaximum(maxRadius);
  h->GetXaxis()->SetTitle( "layer" );
  h->GetYaxis()->SetTitle( "radius" );
  h->Draw();

  h->GetYaxis()->SetTitleOffset( 1.7 );

  tg->SetMarkerStyle(20);
  tg->SetLineColor(1);
  tg->SetMarkerColor(1);
  tg->Draw("P");

  pdfDocument.Add( cv );

  // print
  Stream::PrintVector( "radius", layerRadius, "%.2f" );
  return pdfFile;
}
