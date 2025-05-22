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

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString DeltaZ_qa( TString tag = TString() )
{

  set_style( false );

  if( tag.IsNull() ) tag = "_acts_full_notpc_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/qa*.root", tag.Data() );
  // const TString inputFile = Form( "DST/CONDOR_powerlaw_micromegas/dst_reco%s/qa*.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaZ_qa%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ_qa%s.root", tag.Data() );

  std::cout << "DeltaZ - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto h2d = static_cast<TH2*>(fileManager.GetHistogram( "h_QAG4SimulationDistortions_deltaz_layer" ));

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    if( !is_tpc( idet ) ) continue;
 
    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );
   
    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      const int layerIndex = firstLayer[idet] + ilayer;
      const int binIndex = layerIndex + 1;
      
      const auto hname = Form( "h_%i", layerIndex );
      std::unique_ptr<TH1> h( h2d->ProjectionY( hname, binIndex, binIndex ) );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );

      cv->cd( ilayer+1 );
      h->SetMaximum( 1.2*h->GetMaximum() );
      h->GetXaxis()->SetRangeUser( -0.5, 0.5 );
      h->Draw();

      // fit
      const auto entries( h->GetEntries() );
      std::cout << "DeltaZ - layer: " << layerIndex << " entries: " << entries << std::endl;
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

            auto mean = f->GetParameter(1);
            auto meanError = f->GetParError(1);
            Draw::PutText( 0.2, 0.8, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );

            auto rms = h->GetRMS();
            auto error = f->GetParError(2);
            Draw::PutText( 0.2, 0.7, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

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
    h->SetMaximum(800);
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-cluster) (#mum)" );
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
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-cluster) (#mum)" );
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

    std::cout << "DeltaZ - maximum: " << maximum << std::endl;

  }

  return pdfFile;

}
