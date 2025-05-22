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

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
float DeltaPhi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
double function( double* x, double* )
{
  const double r = x[0];

  static constexpr double rmin_tpc = 20;
  static constexpr double rmax_tpc = 78;
  static constexpr double rlength_tpc = rmax_tpc - rmin_tpc;
  static constexpr double deltarphi_max = 1.5;

  // return -deltarphi_max*std::sin( M_PI*(r-rmin_tpc)/rlength_tpc );
  return -deltarphi_max*std::cos( 2*M_PI*(r-rmin_tpc)/rlength_tpc );

}

//____________________________________________________________________________
TString DeltaRPhi_rphi_distortions( TString tag = TString() )
{
  set_style( false );

  constexpr float max_residual = 5;

  // if( tag.IsNull() ) tag = "_1k_drphi_flat_full_notpc_noouter" ;
  if( tag.IsNull() ) tag = "_1k_drphi_realistic_truth_notpc_noouter" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_distortions%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi_rphi_distortions - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi_rphi_distortions - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  constexpr bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_r*DeltaPhi(_clusters._trk_phi - _clusters._phi)" );
  const TString var2d = Form( "%s:_clusters._trk_r", var.Data() );
  const TString var2dLayer = Form( "%s:_clusters._layer", var.Data() );
  const TCut detector_cut( "_clusters._layer >= 7 && _clusters._layer < 55" );
  const TCut momentum_cut;
  const TCut residual_cut( Form( "fabs(%s)<%f", var.Data(), max_residual ) );

  auto h2 = new TH2F( "h2", "h2", 30, 20, 80, 100, -2.5, 2.5 );
  Utils::TreeToHisto( tree, "h2", var2d, momentum_cut && detector_cut && residual_cut, false );
  h2->SetTitle(TString());
  h2->GetXaxis()->SetTitle( "r (cm)" );
  h2->GetYaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

  auto hlayer = new TH2F( "hlayer", "", nLayers_tpc, firstLayer_tpc, firstLayer_tpc+nLayers_tpc, 100, -max_residual, max_residual );
  Utils::TreeToHisto( tree, hlayer->GetName(), var2dLayer, momentum_cut && detector_cut && residual_cut, false );

  auto cv = new TCanvas( "cv0", "cv0", 800, 800 );
  Draw::DivideCanvas( cv, nLayers_tpc, false );

  // TGraphs
  auto tg = new TGraphErrors();
  auto tg_fit = new TGraphErrors();
  auto tg_fit1 = new TGraphErrors();
  for( int i = 0; i < nLayers_tpc; ++i )
  {

    int layerIndex = firstLayer_tpc + i;

    const auto hname = Form( "h_%i", layerIndex );
    auto h( hlayer->ProjectionY( hname, i+1, i+1 ) );
    h->SetTitle( hname );
    h->SetLineColor( 1 );
    h->SetMarkerColor( 1 );
    h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );
    h->GetXaxis()->SetRangeUser( -max_residual, max_residual );

    cv->cd( i+1 );
    h->Draw();

    if( true )
    {
      double r = radius[layerIndex];
      const double dr = ::function( &r, nullptr );
      tg->SetPoint( i, r, dr );

    } else {

      const double r = radius[layerIndex];
      tg->SetPoint( i, r, 0 );

    }

    // fit
    if( do_fit && h->GetEntries() )
    {
        const auto result = std::min( Fit( h ), Fit_box( h ) );
        auto f = result._function;
        f->Draw("same");
        auto h = f->GetHistogram();
        auto mean = f->GetParameter(1);
        auto error = f->GetParError(1);
        auto rms = h->GetRMS();
        Draw::PutText( 0.2, 0.8, Form( "#LTr#Delta#phi#GT = %.3g #pm %.3g #mum", mean*1e4, rms*1e4 ) );

        tg_fit1->SetPoint( i, radius[layerIndex], mean );
        tg_fit1->SetPointError( i, 0, rms );

        tg_fit->SetPoint( i, radius[layerIndex], mean );
        tg_fit->SetPointError( i, 0, error );
    }

  }

  pdfDocument.Add( cv );

  {
    const TString cvName( "cv" );
    TCanvas* cv( new TCanvas( cvName, cvName, 800, 800 ) );
    cv->SetLeftMargin( 0.14 );
    h2->Draw();

    if( do_fit)
    {
      tg_fit1->SetMarkerStyle(20);
      tg_fit1->SetLineColor(2);
      tg_fit1->SetMarkerColor(2);
      tg_fit1->Draw("P");
    }

    auto f = new TF1( "f", ::function, 20, 80, 0 );
    f->SetLineColor( 4 );
    f->Draw("same" );

    pdfDocument.Add( cv );
  }


  // draw tgraphs
  if( do_fit )
  {
    gStyle->SetOptStat();

    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    tg->SetMarkerStyle( 20 );
    tg->SetMarkerColor( 4 );
    tg->SetLineColor( 4 );
    tg->GetXaxis()->SetTitle( "r (cm)" );
    tg->GetYaxis()->SetTitle( "#Deltar (cm)" );
    tg->SetMinimum( -2.5 );
    tg->SetMaximum( 2.5 );
    tg->Draw("AP");

    tg_fit->SetMarkerStyle( 20 );
    tg_fit->SetMarkerColor( 2 );
    tg_fit->SetLineColor( 2 );
    tg_fit->Draw("P");
    cv->Update();

    auto legend = new TLegend( 0.16, 0.82, 0.60, 0.93, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();
    legend->AddEntry( tg, "input distortion", "AP" );
    legend->AddEntry( tg_fit, "reconstructed distortion", "AP" );

    pdfDocument.Add( cv );
  }

  if( do_fit )
  {

    auto pulls = new TH1F( "pulls", "pulls", 100, -5, 5 );
    pulls->GetXaxis()->SetTitle( "pulls" );
    for( int i = 0; i < tg->GetN(); ++i )
    {
      double x, y;
      tg->GetPoint( i, x, y );

      double x_fit, y_fit;
      tg_fit->GetPoint( i, x_fit, y_fit );

      double error = tg_fit->GetErrorY( i );
      pulls->Fill( (y_fit - y)/error );
    }

    TCanvas* cv( new TCanvas( "cv2", "cv2", 800, 800 ) );
    pulls->Draw();
    pdfDocument.Add( cv );

  }

  return pdfFile;

}
