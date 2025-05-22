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
#include "Fit.C"

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhi_truth_new( TString tag = TString() )
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.003, 0.01, 0.2, 0.2, 0.2, 0.5, 0.5};

  if( tag.IsNull() ) tag = "_genfit_truth_realistic-new";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_truth%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhi_truth%s.root", tag.Data() );

  std::cout << "DeltaRPhi_truth - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "m_tracks.clusters[].trk_r*delta_phi(m_tracks.clusters[].trk_phi - m_tracks.clusters[].truth_phi)" );
  const TString var2d = Form( "%s:m_tracks.clusters[].layer", var.Data() );
  const TCut momentum_cut = ("m_tracks.pt>0.5" );
  const TCut pattern_cut(
    "m_tracks.truth_pt>0.5"
    "&&m_tracks.nclusters_mvtx>0"
    "&&m_tracks.nclusters_intt>=2"
    // "&&m_tracks.nclusters_micromegas==2"
    );
  
  //   const TCut pattern_cut("m_tracks.nclusters_micromegas==0");
    
  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );

  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    // skip detector 6, which is z segmented micromegas
    if( idet == 6 ) continue;

    const TString hname( Form( "deltarphi_%i_0", idet ) );
    const TCut layer_cut( Form( "m_tracks.clusters[].layer==%i", firstLayer[idet] ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut&&pattern_cut, false );
      max_det_residual[idet] = 8*h1->GetRMS() + std::abs( h1->GetMean() );;
    }
  }

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    // skip detector 6, which is z segmented micromegas
    if( idet == 6 ) continue;

    const TCut detector_cut( Form( "m_tracks.clusters[].layer>=%i &&m_tracks.clusters[].layer<%i ", firstLayer[idet], firstLayer[idet] + nLayers[idet] ) );

    const TString hname( Form( "deltarphi_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut&&detector_cut, false );

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
      h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-truth} (cm)" );
      h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );

      cv->cd( ilayer+1 );
      h->SetMaximum( 1.2*h->GetMaximum() );
      h->Draw();

      // fit
      const auto entries( h->GetEntries() );
      std::cout << "DeltaRPhi_truth - layer: " << layerIndex << " entries: " << entries << std::endl;
      if( entries )
      {
        if( do_fit )
        {
          const auto result = std::min( Fit( h.get() ), Fit_box( h.get() ) );
          if( result._valid )
          {
            auto f = result._function;
            f->Draw("same");
            auto h = f->GetHistogram();

            Draw::PutText( 0.2, 0.85, Form( "entries = %.0f", entries ) );            
            
            auto mean = f->GetParameter(1);
            auto meanError = f->GetParError(1);
            Draw::PutText( 0.2, 0.75, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );

            auto rms = h->GetRMS();
            auto error = f->GetParError(2);
            Draw::PutText( 0.2, 0.65, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

            tgl->SetPoint( layerIndex, layerIndex, rms*1e4 );
            tgl->SetPointError( layerIndex, 0, error*1e4 );

            tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
            tg->SetPointError( layerIndex, 0, error*1e4 );
          } else {
            std::cout << "DeltaRPhi_truth - skipping layer " << layerIndex << " (failed fit)" << std::endl;
          }
        } else {
          auto mean = h->GetMean();
          auto meanError = h->GetMeanError();
          Draw::PutText( 0.2, 0.8, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );

          auto rms = h->GetRMS();
          auto error = h->GetRMSError();
          Draw::PutText( 0.2, 0.7, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

          tgl->SetPoint( layerIndex, layerIndex, rms*1e4 );
          tgl->SetPointError( layerIndex, 0, error*1e4 );

          tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
          tg->SetPointError( layerIndex, 0, error*1e4 );
        }
      }

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

      // save in array
      h_array[layerIndex] = std::move(h);
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
    h->SetMaximum(1.2*Utils::GetMaximum(tgl));
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
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
    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 20, 80 ) );
    h->SetMinimum(0);
    // h->SetMaximum(160);
    h->SetMaximum(1.2*Utils::GetMaximum(tg));
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
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
  for( auto&& h:h_array) { if(h) h->Write(); }
  tgl->Write();
  tg->Write();
  output->Close();

  if( true )
  {

    // get the maximum in the tpc
    double maximum = 0;
    for( int ip = 0; ip < tg->GetN(); ++ip )
    {
      double r, sigma;
      tg->GetPoint( ip, r, sigma );
      if( r < rmin_tpc || r > rmax_tpc ) continue;
      if( sigma > maximum ) maximum = sigma;
    }

    std::cout << "DeltaRPhi - maximum: " << maximum << std::endl;

  }

  return pdfFile;

}
