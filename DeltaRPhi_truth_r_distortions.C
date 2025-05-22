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
double function( double* x, double* par )
{
  const double ta = x[0];
  return par[0]+ta*par[1];
}

//____________________________________________________________________________
TString DeltaRPhi_truth_r_distortions( TString tag = TString() )
{

  set_style( false );
  gStyle->SetPadLeftMargin(0.13);

  constexpr float max_residual = 4;
  constexpr float max_alpha = 0.5;

//   if( tag.IsNull() ) tag = "_1k_dr_realistic_truth_notpc_ross" ;
//   if( tag.IsNull() ) tag = "_1k_dr_realistic_truth_notpc" ;
//   if( tag.IsNull() ) tag = "_1k_realistic_full_notpc_noouter" ;
//   if( tag.IsNull() ) tag = "_1k_dr_realistic_truth_notpc_nphi1k" ;
//   if( tag.IsNull() ) tag = "_1k_dr_realistic_truth_notpc_single_nphi1k" ;
  if( tag.IsNull() ) tag = "_1k_dr_realistic_truth_notpc_noouter" ;
//   if( tag.IsNull() ) tag = "_1k_realistic_truth_notpc_noouter" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

//   if( tag.IsNull() ) tag = "_dr_realistic_truth_notpc_noouter" ;
//   const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_truth_r_distortions%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi2D - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi2D - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  constexpr bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_r*DeltaPhi(_clusters._trk_phi - _clusters._phi)" );
  const TString var2d = Form( "%s:tan( _clusters._truth_alpha )", var.Data() );
  const TCut momentum_cut;
  const TCut alpha_cut( Form( "fabs( _clusters._trk_alpha ) < %f", max_alpha ) );
  const TCut residual_cut( Form( "fabs(%s)<%f", var.Data(), max_residual ) );

//   constexpr int nslices = 6;
//   constexpr std::array<int, nslices+1> layers = { 7, 15, 23, 31, 39, 47, 55 };

  constexpr int nslices = 12;
  constexpr std::array<int, nslices+1> layers = { 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55 };

  TCanvas* cv( new TCanvas( "cv", "cv", 800, 800 ) );
  //   Draw::DivideCanvas( cv, nslices, false );
  cv->Divide( 4, 3, 0, 0 );

  // TGraphs
  auto tg = new TGraphErrors();
  auto tg_fit = new TGraphErrors();
  for( int i = 0; i < nslices; ++i )
  {
    const TCut detector_cut( Form( "_clusters._layer >= %i && _clusters._layer < %i", layers[i], layers[i+1] ) );

    auto hname = Form( "h%i", i );
    auto h = new TH2F( hname, hname, 100, -max_alpha, max_alpha, 100, -max_residual, max_residual );
    Utils::TreeToHisto( tree, hname, var2d, momentum_cut && detector_cut && alpha_cut && residual_cut, false );
    h->SetTitle(TString());
    h->GetXaxis()->SetTitle( "tan(#alpha)" );
    h->GetYaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

    auto pname = Form( "p%i", i );
    auto p = new TProfile( pname, pname, 100, -max_alpha, max_alpha );
    Utils::TreeToHisto( tree, pname, var2d, momentum_cut && detector_cut && alpha_cut && residual_cut, false );

    cv->cd( i+1 );
    h->Draw();
    p->SetLineColor(2);
    p->Draw( "same" );

    // calculate mean radius for this slice
    const float r = 0.5*(radius[layers[i]] + radius[layers[i+1]-1]);

    if( false )
    {
      auto fname = Form( "f%i", i );
      auto f = new TF1( fname, ::function, -max_alpha, max_alpha, 2 );
      f->SetParameter( 0, 0 );

      // define the r transformation (cm)
      static constexpr float deltar_max = 1.5;
      const float dr = deltar_max*std::cos( 2*M_PI*(r-rmin_tpc)/rlength_tpc );

      std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta r: " << dr << std::endl;
      f->SetParameter( 1, dr );
      f->SetLineColor(4);
      f->Draw("same");

      tg->SetPoint( i, r, dr );

    } else {

      tg->SetPoint( i, r, 0 );

    }

    if( do_fit && h->GetEntries() )
    {
      auto fname = Form( "f_fit%i", i );
      auto f = new TF1( fname, ::function, -0.3, 0.3, 2 );
      f->SetParameter( 0, 0 );
      f->SetParameter( 1, 0 );

      h->Fit( f, "0Q" );
      f->SetLineColor( 2 );
      f->Draw("same");

      tg_fit->SetPoint( i, r, f->GetParameter(1) );
      tg_fit->SetPointError( i, 0, f->GetParError(1) );

    }

  }

  pdfDocument.Add( cv );

  // draw tgraphs
  if( do_fit )
  {
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

  return pdfFile;

}
