#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString DeltaZ( TString tag = TString() )
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // const std::array<Float_t, nDetectors> maxDetResidual = { 0.005, 0.015, 0.05, 0.05, 0.05, 0.05 };
  // const std::array<Float_t, nDetectors> maxDetResidual = { 0.003, 0.01, 0.09, 0.09, 0.09, 0.03 };
  // const std::array<Float_t, nDetectors> maxDetResidual = { 0.003, 0.01, 0.15, 0.15, 0.15, 0.15 };
  const std::array<Float_t, nDetectors> maxDetResidual = { 0.003, 1, 0.3, 0.3, 0.3, 0.15 };

  auto maxResidual = *std::max_element( maxDetResidual.cbegin(), maxDetResidual.cend() )/3;

  // pdf output
  if( tag.IsNull() ) tag = "_5k_full_notpc_noouter" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaZ%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ%s.root", tag.Data() );

  std::cout << "DeltaZ - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const Bool_t doFit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_z - _clusters._z" );
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // save all histograms
  std::array<TH1*, nLayersTotal> h_array;

  // loop over detectors
  for( Int_t idet = 0; idet < nDetectors; ++idet )
  {

    const TString hName( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hName, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -maxDetResidual[idet], maxDetResidual[idet] ) );
    Utils::TreeToHisto( tree, hName, var2d, TCut(), kFALSE );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    Draw::DivideCanvas( cv, nLayers[idet], kFALSE );

    // loop over layers
    for( Int_t ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      Int_t layerIndex = firstLayer[idet] + ilayer;
      const auto hName = Form( "h_%i", layerIndex );
      TH1* h = h2d->ProjectionY( hName, ilayer+1, ilayer+1 );
      h->SetTitle( hName );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{trk-clus} (cm)" );
      h->GetXaxis()->SetRangeUser( -maxDetResidual[idet], maxDetResidual[idet] );
      h->SetMaximum( h->GetMaximum()*1.2 );

      cv->cd( ilayer+1 );
      h->Draw();

      // fit
      if( doFit )
      {
        auto f = Fit( h );
        Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", f->GetParameter(2)*1e4, f->GetParError(2)*1e4 ) );
        tg->SetPoint( layerIndex, radius[layerIndex], f->GetParameter(2)*1e4 );
        tg->SetPointError( layerIndex, 0, f->GetParError(2)*1e4 );

      } else {

        Draw::PutText( 0.2, 0.8, Form( "RMS = %.3g #pm %.3g #mum", h->GetRMS()*1e4, h->GetRMSError()*1e4 ) );
        tg->SetPoint( layerIndex, radius[layerIndex], h->GetRMS()*1e4 );
        tg->SetPointError( layerIndex, 0, h->GetRMSError()*1e4 );

      }

      // save in array
      h_array[layerIndex] = h;

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv );

  }

  // TGraph
  auto cv = new TCanvas( "cvtg", "cvtg", 800, 800 );
  cv->SetLeftMargin( 0.15 );

  auto h = new TH1F( "dummy", "", 100, 0, 90 );
  h->SetMinimum(0);
  h->SetMaximum(maxResidual*1e4);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-clus) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.7 );
  h->Draw();

  tg->SetMarkerStyle(20);
  tg->SetLineColor(1);
  tg->SetMarkerColor(1);
  tg->Draw("P");

  pdfDocument.Add( cv );

  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:h_array) { h->Write(); }
  tg->Write();
  output->Close();

  return pdfFile;

}
