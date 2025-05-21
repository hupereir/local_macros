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
TString DeltaZ()
{

  set_style( false );

  std::array<float, nDetectors> max_det_residual = { 0.003, 1, 0.5, 0.5, 0.5, 0.1, 0.1 };
  // std::array<float, nDetectors> max_det_residual = { 0.003, 1, 2, 2, 2, 0.1, 0.1 };
  
  // input files
  // input files
  // const TString tag = "_acts_full_notpc_nodistortion-new";
  const TString tag = "_genfit_full_notpc_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco*_?.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/DeltaZ%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ%s.root", tag.Data() );

  std::cout << "DeltaZ - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._trk_z - _tracks._clusters._z" );
  const TString var2d = Form( "%s:_tracks._clusters._layer", var.Data() );
  const TCut momentum_cut("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    );
//   const TCut momentum_cut;
//   const TCut pattern_cut;

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );
  
  // create TGraph to store resolution vs layer
  auto tgm = new TGraphErrors();
  tgm->SetName( "residuals mean" );

  // create TGraph to store resolution vs radius
  auto tgml = new TGraphErrors();
  tgml->SetName( "residuals_layers mean" );
  
  // optimize max residual
  if( true )
  {
    for( int idet = 0; idet < nDetectors; ++idet )
    {
      // skip detector 5, which is phi segmented micromegas
      if( idet == 5 ) continue;
      
      const TString hname( Form( "deltaz_%i_0", idet ) );
      const TCut layer_cut( Form( "_tracks._clusters._layer==%i", firstLayer[idet] ) );
      
      for( int i=0; i<3; ++i )
      {
        std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
        Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut&&pattern_cut, false );
        max_det_residual[idet] = 5*h1->GetRMS() + std::abs( h1->GetMean() );;
      }
    }
  }
  
  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    // skip detector 5, which is phi segmented micromegas
    if( idet == 5 ) continue;

    const TCut detector_cut( Form( "_tracks._clusters[]._layer>=%i &&_tracks._clusters[]._layer<%i ", firstLayer[idet], firstLayer[idet] + nLayers[idet] ) );

    const TString hname( Form( "deltaz_%i", idet ) );
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
      const auto entries( h->GetEntries() );
      std::cout << "DeltaRPhi_truth - layer: " << layerIndex << " entries: " << entries << std::endl;
      if( !entries ) continue;

      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{trk-clus} (cm)" );
      h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
      h->SetMaximum( h->GetMaximum()*1.2 );

      cv->cd( ilayer+1 );
      h->Draw();

      // fit
      if( do_fit )
      {
        const auto result = std::min( { Fit( h.get() ), Fit_box( h.get() ) } );
        if( result._valid )
        {
          auto f = result._function;
          f->Draw("same");
          auto h = f->GetHistogram();
          
          Draw::PutText( 0.2, 0.85, Form( "entries = %.0f", entries ) );            
          
          auto mean = f->GetParameter(1);
          auto meanError = f->GetParError(1);
          Draw::PutText( 0.2, 0.75, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );
          
          tgml->SetPoint( layerIndex, layerIndex, mean );
          tgml->SetPointError( layerIndex, 0, meanError );
          
          tgm->SetPoint( layerIndex, radius[layerIndex], mean );
          tgm->SetPointError( layerIndex, 0, meanError );
          
          auto rms = h->GetRMS();
          auto error = f->GetParError(2);
          Draw::PutText( 0.2, 0.65, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );
          
          tgl->SetPoint( layerIndex, layerIndex, rms*1e4 );
          tgl->SetPointError( layerIndex, 0, error*1e4 );
          
          tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
          tg->SetPointError( layerIndex, 0, error*1e4 );
        } else {
          std::cout << "DeltaZ - skipping layer " << layerIndex << " (failed fit)" << std::endl;
        }
        
      } else {
        auto mean = h->GetMean();
        auto meanError = h->GetMeanError();
        Draw::PutText( 0.2, 0.8, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );
        
        tgml->SetPoint( layerIndex, layerIndex, mean );
        tgml->SetPointError( layerIndex, 0, meanError );

        tgm->SetPoint( layerIndex, radius[layerIndex], mean );
        tgm->SetPointError( layerIndex, 0, meanError );

        auto rms = h->GetRMS();
        auto error = h->GetRMSError();
        Draw::PutText( 0.2, 0.7, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );
        
        tgl->SetPoint( layerIndex, layerIndex, rms*1e4 );
        tgl->SetPointError( layerIndex, 0, error*1e4 );
        
        tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
        tg->SetPointError( layerIndex, 0, error*1e4 );
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

  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum(tgl));
    h->SetMaximum(900);
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-cluster) (#mum)" );
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
    h->SetMaximum(1.2*Utils::GetMaximum(tg));
    h->SetMaximum(900);
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#sigma_{#Deltaz} (track-cluster) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tg->SetMarkerStyle(20);
    tg->SetLineColor(1);
    tg->SetMarkerColor(1);
    tg->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // static constexpr float max_offset = 1.;
  static constexpr float max_offset = 0.05;

  // mean vs layer
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(-max_offset);
    h->SetMaximum(max_offset);
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#LTr.#Deltaz#GT (track-cluster) (cm)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tgml->SetMarkerStyle(20);
    tgml->SetLineColor(1);
    tgml->SetMarkerColor(1);
    tgml->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // Mean vs radius
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(-max_offset);
    h->SetMaximum(max_offset);
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#LTr.#Deltaz#GT (track-cluster) (cm)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tgm->SetMarkerStyle(20);
    tgm->SetLineColor(1);
    tgm->SetMarkerColor(1);
    tgm->Draw("P");

    pdfDocument.Add( cv.get() );
  }

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

    std::cout << "DeltaZ - maximum: " << maximum << std::endl;

  }
  
  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:h_array) { if(h) h->Write(); }
  tgl->Write();
  tg->Write();
  tgml->Write();
  tgm->Write();
  output->Close();

  return pdfFile;

}
