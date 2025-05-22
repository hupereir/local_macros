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
TString DeltaZ_Tony()
{

  set_style( false );

  // const TString tag = "_corrected";
  // const TString tag = "_corrected_notpc";
  const TString tag = "_corrected_notpc-rotated9";

  // run list
  const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };

//   // run list
//   const std::vector<int> runlist = {49714, 49715, 49716, 49717, 49718, 49719, 49720, 49721, 49722 };
//   const TString postfix = "_6X6";

//   // run list
//   const std::vector<int> runlist = {50015, 50016, 50017, 50018, 50019, 50020};
//   const TString postfix = "_28x28";

//   // run list
//   const std::vector<int> runlist = { 50450, 50451, 50452, 50454 };
//   const TString postfix = "_56x56";  // input files
  const TString tag = "_corrected_notpc";
  const TString inputFile = Form( "DST/CONDOR_CombinedDataReconstruction%s/TrackResiduals-*-full.root", tag.Data() );

//   const TCut charge_cut;
//   const TString postfix;

//   const TCut charge_cut("charge<0");
//   const TString postfix = "_neg";

  const TCut charge_cut("charge>0");
  const TString postfix = "_pos";


  const TString pdfFile = Form( "Figures/DeltaZ_Tony%s%s.pdf", tag.Data(), postfix.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ_Tony%s%s.root", tag.Data(), postfix.Data());

  std::cout << "DeltaZ_Tony - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_Tony - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  //   const double max_dz = 2;
  const double max_dz = 5;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "residualtree" );
  if( !tree ) return TString();

  // variable names
  const TString var( "stategz - clusgz" );
  const TString var2d = Form( "%s:cluslayer", var.Data() );
  const TCut track_cut(
    "m_pt>0.2 && quality<100"
    "&& m_ntpc>20 && m_nmaps>2 && m_nintt>1"
    "&& nmms>0");

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
      const TCut layer_cut( Form( "cluslayer==%i", firstLayer[idet] ) );

      for( int i=0; i<3; ++i )
      {
        std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
        Utils::TreeToHisto( tree, hname, var, track_cut&&charge_cut&&layer_cut, false );
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

    const TCut layer_cut( Form( "cluslayer>=%i &&cluslayer<%i ", firstLayer[idet], firstLayer[idet] + nLayers[idet] ) );

    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var2d, track_cut&&charge_cut&&layer_cut, false );

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
    // rms vs layer
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

  // static constexpr float max_offset = 1.;
  static constexpr float max_offset = 0.05;

  // mean vs layer
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    {
      // rphi residuals vs layer
      const TString var = "stategz-clusgz";
      const TString var2d = Form( "%s:cluslayer", var.Data() );

      const auto hname = "h_drphi_2d";
      auto h = new TH2F( hname, "", 60, 0, 60, 100, -max_dz, max_dz );
      Utils::TreeToHisto( tree, hname, var2d, track_cut&&charge_cut, false );

      h->GetXaxis()->SetTitle( "layer" );
      h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );
      h->SetStats(0);
      h->Draw( "colz" );
    }

    tgml->SetMarkerStyle(20);
    tgml->SetLineColor(1);
    tgml->SetMarkerColor(1);
    tgml->Draw("P");

    gPad->Update();
    Draw::HorizontalLine(gPad, 0)->Draw();

    pdfDocument.Add( cv.get() );
  }

  // mean vs radius
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgr", "cvtgr", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    {
      // rphi residuals vs layer
      const TString var = "stategz-clusgz";
      const TString var2d = Form( "%s:fabs(clusgr)", var.Data() );

      const auto hname = "h_drphi_r_2d";
      auto h = new TH2F( hname, "", 100, 0, 90, 100, -max_dz, max_dz );
      Utils::TreeToHisto( tree, hname, var2d, track_cut&&charge_cut, false );

      h->GetXaxis()->SetTitle( "r_{clus} (cm)" );
      h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );
      h->SetStats(0);
      h->Draw( "colz" );
    }

    tgm->SetMarkerStyle(20);
    tgm->SetLineColor(1);
    tgm->SetMarkerColor(1);
    tgm->Draw("P");

    gPad->Update();
    Draw::HorizontalLine(gPad, 0)->Draw();

    pdfDocument.Add( cv.get() );
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
