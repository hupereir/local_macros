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
TString PullsZ_truth( TString tag = TString() )
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::array<Float_t, nDetectors> max_det_pull = { 5, 5, 5, 5, 5 };

  // input files
  if( tag.IsNull() ) tag = "_2k_realistic_full_nominal_old" ;
  const TString inputFile = Form( "DST/dst_eval%s*.root", tag.Data() );

  const TString pdfFile = Form( "Figures/PullsZ_truth%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/PullsZ_truth%s.root", tag.Data() );

  std::cout << "PullsZ_truth - inputFile: " << inputFile << std::endl;
  std::cout << "PullsZ_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "(_tracks._clusters[]._trk_z - _tracks._clusters[]._truth_z)/_tracks._clusters[]._trk_z_error" );
  const TString var2d = Form( "%s:_tracks._clusters[]._layer", var.Data() );
  const TCut momentum_cut;

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "pulls" );

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "pulls_layers" );

  // optimize max pull
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    const TString hname( Form( "deltaz_%i_0", idet ) );
    const TCut layer_cut( Form( "_tracks._clusters[]._layer==%i", firstLayer[idet]+1 ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_pull[idet], max_det_pull[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut, false );
      max_det_pull[idet] = 8*h1->GetRMS();
    }

  }

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_pull[idet], max_det_pull[idet] ) );
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
      std::unique_ptr<TH1> h( h2d->ProjectionY( hname, ilayer+1, ilayer+1 ) );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{track-truth}/#sigmaz" );
      h->GetXaxis()->SetRangeUser( -max_det_pull[idet], max_det_pull[idet] );

      cv->cd( ilayer+1 );
      h->SetMaximum( 1.2*h->GetMaximum() );
      h->Draw();

      // fit
      const auto entries( h->GetEntries() );
      if( entries )
      {
        if( do_fit )
        {
          auto result = std::min( { Fit( h.get() ), Fit_box( h.get() ) } );
          if( result._valid )
          {
            auto f = result._function;
            f->Draw("same");
            auto h = f->GetHistogram();
            auto rms = h->GetRMS();
            auto error = f->GetParError(2);

            Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g", rms, error ) );
            Draw::PutText( 0.2, 0.7, Form( "entries: %.0f", entries ) );

            tgl->SetPoint( layerIndex, layerIndex, rms );
            tgl->SetPointError( layerIndex, 0, error );

            tg->SetPoint( layerIndex, radius[layerIndex], rms );
            tg->SetPointError( layerIndex, 0, error );

          } else {

            std::cout << "PullsZ_truth - skipping layer " << layerIndex << " (failed fit)" << std::endl;

          }

        } else {

          auto rms = h->GetRMS();
          auto error = h->GetRMSError(2);
          Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g", rms, error ) );

          tgl->SetPoint( layerIndex, layerIndex, rms );
          tgl->SetPointError( layerIndex, 0, error );

          tg->SetPoint( layerIndex, radius[layerIndex], rms );
          tg->SetPointError( layerIndex, 0, error );

        }
      }

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

      // save in array
      h_array[layerIndex] = std::move( h );

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
    h->SetMaximum( 1.2*Utils::GetMaximum( tgl ) );
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "pull_{#Deltaz} (track-truth)" );
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

    // std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(0);
    h->SetMaximum( 1.2*Utils::GetMaximum( tg ) );
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "pull_{#Deltaz} (track-truth)" );
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
