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

#include <trackbase/TrkrClusterContainer.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

//____________________________________________________________________________
TString make_run_label( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "run %i", runlist[0] );
  return Form( "runs %i-%i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
TString make_run_postfix( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "_%08i", runlist[0] );
  return Form( "_%08i-%08i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString Phi_raw()
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 3.14,3.14,3.14,3.14,3.14,3.14,3.14};

  // input files
  const TString tag = "_CombinedDataReconstruction_notpc";
  TString run_label;
  TString postfix;

  FileManager fileManager;

  {
    const std::vector<int> runlist = { 53534 };
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );
  }

  const TString pdfFile = Form( "Figures/Phi_raw%s%s.pdf", tag.Data(), postfix.Data() );
  const TString rootFile  = Form( "Rootfiles/Phi_raw%s%s.root", tag.Data(), postfix.Data() );

  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );


  // file manager
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  // entries
  std::cout << "DeltaRPhi_raw - entries: " << tree->GetEntries() << std::endl;
  Utils::max_entries = 10000;

  // variable names
  const TString var( "_tracks._clusters._trk_phi" );
  const TString var2d = Form( "%s:_tracks._clusters._layer", var.Data() );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._pt>0.5"
    "&&_tracks._nclusters_tpc>20"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._clusters._trk_r > 0"
    );

//   const TCut momentum_cut;
//   const TCut pattern_cut;

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName("residuals");

  // create TGraph to store resolution vs layer
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );

  // create TGraph to store mean vs layer
  auto tgm = new TGraphErrors();
  tgm->SetName( "residuals mean" );

  // create TGraph to store mean vs radius
  auto tgml = new TGraphErrors();
  tgml->SetName( "residuals_layers mean" );

  // save all histograms
  std::array<std::unique_ptr<TH1>, nLayersTotal> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // skip detector 6, which is z segmented micromegas
    if( idet == 6 ) continue;

    const TCut layer_cut( Form( "_tracks._clusters[]._layer>=%i &&_tracks._clusters[]._layer<%i ", firstLayer[idet], firstLayer[idet] + nLayers[idet] ) );

    const TString hname( Form( "deltarphi_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut&&layer_cut, false );

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
      std::cout << "DeltaRPhi - layer: " << layerIndex << " entries: " << entries << std::endl;
      if( !entries ) continue;

      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );
      h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );

      cv->cd( ilayer+1 );
      h->SetMaximum( 1.2*h->GetMaximum() );
      h->Draw();

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

  // RMS vs layer
  const double max_rms = 1.2*Utils::GetMaximum(tg);
  // const double max_rms = 350;
  // const double max_rms = 12000;
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(0);
    h->SetMaximum(max_rms);
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

  // RMS vs radius
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    // std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 20, 80 ) );
    h->SetMinimum(0);
    h->SetMaximum(max_rms);
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
    h->GetYaxis()->SetTitle( "#LTr.#Delta#phi#GT (track-cluster) (cm)" );
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
    h->GetYaxis()->SetTitle( "#LTr.#Delta#phi#GT (track-cluster) (cm)" );
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

    std::cout << "DeltaRPhi - maximum: " << maximum << std::endl;

  }

  // save everything in rootfiles
  RootFile output( rootFile );
  for( auto&& h:h_array) { output.Add( h.release() ); }
  output.Add( tgl );
  output.Add( tg );
  output.Add( tgml );
  output.Add( tgm );

  return pdfFile;

}
