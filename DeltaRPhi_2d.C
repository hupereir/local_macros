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
  if( phi >= 2*M_PI ) return phi - 2*M_PI;
  else if( phi <= -2*M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
double function( double* x, double* )
{
  const double r = x[0];

  static constexpr double rmin_tpc = 20;
  static constexpr double rmax_tpc = 78;
  static constexpr double rlength_tpc = rmax_tpc - rmin_tpc;
  static constexpr double deltarphi_max = 0.3;

  return -deltarphi_max*std::sin( M_PI*(r-rmin_tpc)/rlength_tpc );

}

//____________________________________________________________________________
TString DeltaRPhi_2d( TString tag = TString() )
{

  set_style( false );


  constexpr float max_residual = 1.2;

  if( tag.IsNull() ) tag = "_1k_drphi_notpc_noouter_flat" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_2d%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi2D - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi2D - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_r*DeltaPhi(_clusters._trk_phi - _clusters._phi)" );
  const TString var2d = Form( "%s:_clusters._trk_r", var.Data() );
  const TCut detector_cut( "_clusters._layer >= 7 && _clusters._layer < 55" );
  const TCut momentum_cut;
  const TCut residual_cut( Form( "fabs(%s)<%f", var.Data(), max_residual ) );

  auto h = Utils::TreeToHisto( tree, "h", var2d, momentum_cut && detector_cut && residual_cut, true );
  h->SetTitle(TString());
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

  const TString cvName( "cv" );
  TCanvas* cv( new TCanvas( cvName, cvName, 800, 800 ) );
  h->Draw();

  auto f = new TF1( "f", ::function, 20, 80, 0 );
  f->SetLineColor( 2 );
  f->Draw("same" );

  pdfDocument.Add( cv );

  return pdfFile;

}
