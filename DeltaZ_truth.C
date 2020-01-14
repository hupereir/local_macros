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
TString DeltaZ_truth( TString tag = TString() )
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::array<Float_t, nDetectors> max_det_residual = { 0.005, 0.015, 0.25, 0.25, 0.4, 0.4 };

//   if( tag.IsNull() ) tag = "_1k_truth_notpc_noouter" ;
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  if( tag.IsNull() ) tag = "_realistic_truth_notpc_single_nz500" ;
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaZ_truth%s_highpt.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ_truth%s_highpt.root", tag.Data() );

  std::cout << "DeltaZ_truth - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool doFit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._trk_z - _clusters._truth_z" );
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );
  TCut momentum_cut( "_clusters._truth_pt>6" );
  // TCut momentum_cut;

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );

  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    const TString hname( Form( "deltaz_%i_0", idet ) );
    const TCut layer_cut( Form( "_clusters._layer==%i", firstLayer[idet]+1 ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut, false );
      max_det_residual[idet] = 5*h1->GetRMS();
    }

    max_det_residual[idet] = 8*max_det_residual[idet]/5;

  }

  const auto max_residual = *std::max_element( max_det_residual.cbegin(), max_det_residual.cend() )/5;

  // save all histograms
  std::array<TH1*, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 50, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var2d, momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      auto h( h2d->ProjectionY( hname, ilayer+1, ilayer+1 ) );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{track-truth} (cm)" );
      h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );

      cv->cd( ilayer+1 );
      h->Draw();

      // fit
      if( h->GetEntries() )
      {
        if( doFit )
        {
          auto results = { Fit( h ), Fit_box( h ) };
          const auto result = std::min( results );
          auto f = result._function;
          f->Draw("same");
          auto h = f->GetHistogram();
          auto rms = h->GetRMS();
          auto error = f->GetParError(2);

          Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

          tgl->SetPoint( layerIndex, layerIndex, rms*1e4 );
          tgl->SetPointError( layerIndex, 0, error*1e4 );

          tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
          tg->SetPointError( layerIndex, 0, error*1e4 );

        } else {

          tgl->SetPoint( layerIndex, layerIndex, h->GetRMS()*1e4 );
          tgl->SetPointError( layerIndex, 0, h->GetRMSError()*1e4 );

          tg->SetPoint( layerIndex, radius[layerIndex], h->GetRMS()*1e4 );
          tg->SetPointError( layerIndex, 0, h->GetRMSError()*1e4 );

        }
      }

      // save in array
      h_array[layerIndex] = h;

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  // TGraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(0);
    h->SetMaximum(max_residual*1e4);
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tgl->SetMarkerStyle(20);
    tgl->SetLineColor(1);
    tgl->SetMarkerColor(1);
    tgl->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // TGraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(0);
    h->SetMaximum(max_residual*1e4);
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tg->SetMarkerStyle(20);
    tg->SetLineColor(1);
    tg->SetMarkerColor(1);
    tg->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:h_array) { h->Write(); }
  tgl->Write();
  tg->Write();
  output->Close();

  return pdfFile;

}
