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
Float_t DeltaPhi( Float_t phi )
{
  if( phi >= 2*M_PI ) return phi - 2*M_PI;
  else if( phi <= -2*M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
void DeltaRPhi()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // const std::array<Float_t, nDetectors> maxDetResidual = { 0.005, 0.02, 0.1, 0.1, 0.1, 0.1 };
  const std::array<Float_t, nDetectors> maxDetResidual = { 0.005, 0.02, 0.3, 0.3, 0.3, 0.1 };

  // pdf output
  const TString tag = "_5k_full" ;
  // const TString tag = "_5k_full_notpc" ;
  // const TString tag = "_5k_full_notpc_noouter" ;
  // const TString tag = "_full_notpc" ;
  // const TString tag = "_full_notpc_noouter" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaRPhi%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const Bool_t doFit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_trk_r*DeltaPhi(_trk_phi - _phi)" );
  const TString var2d = Form( "%s:_layer", var.Data() );

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();

  // loop over detectors
  for( Int_t idet = 0; idet < nDetectors; ++idet )
  {
    const TString hName( Form( "deltarphi_%i", idet ) );
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
      h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-clus} (cm)" );
      h->GetXaxis()->SetRangeUser( -maxDetResidual[idet], maxDetResidual[idet] );

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

        tg->SetPoint( layerIndex, radius[layerIndex], h->GetRMS()*1e4 );
        tg->SetPointError( layerIndex, 0, h->GetRMSError()*1e4 );

      }

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

  // auto h = new TH1F( "dummy", "", 100, 0, nLayersTotal );
  auto h = new TH1F( "dummy", "", 100, 0, 90 );
  h->SetMinimum(0);
  h->SetMaximum(5e3);
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-clus) (#mum)" );
  h->Draw();

  h->GetYaxis()->SetTitleOffset( 1.7 );

  tg->SetMarkerStyle(20);
  tg->SetLineColor(1);
  tg->SetMarkerColor(1);
  tg->Draw("P");


  pdfDocument.Add( cv );

}
