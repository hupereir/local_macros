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
double function( double* x, double* par )
{
  const double ta = x[0];
  return par[0]+ta*par[1];
}

//____________________________________________________________________________
TString DeltaZ_rdistortions( TString tag = TString() )
{

  set_style( false );


  constexpr float max_residual = 4;
  constexpr float max_beta = 0.5;

  // if( tag.IsNull() ) tag = "_1k_dr_notpc_realistic_ross" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaZ_rdistortions%s.pdf", tag.Data() );

  std::cout << "DeltaZ_rdistortions - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_rdistortions - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_z - _clusters._z" );
  const TString var2d = Form( "%s:_clusters._trk_beta", var.Data() );
  const TCut momentum_cut;
  const TCut beta_cut( Form( "fabs( _clusters._trk_beta ) < %f", max_beta ) );
  const TCut residual_cut( Form( "fabs(%s)<%f", var.Data(), max_residual ) );

//   constexpr int nslices = 6;
//   constexpr std::array<int, nslices+1> layers = { 7, 15, 23, 31, 39, 47, 55 };

  constexpr int nslices = 12;
  constexpr std::array<int, nslices+1> layers = { 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55 };

  const TString cvName( "cv" );
  TCanvas* cv( new TCanvas( cvName, cvName, 800, 800 ) );
  Draw::DivideCanvas( cv, nslices, false );

  for( int i = 0; i < nslices; ++i )
  {
    auto hname = Form( "h%i", i );
    auto h = new TH2F( hname, hname, 100, -max_beta, max_beta, 100, -max_residual, max_residual );

    const TCut detector_cut( Form( "_clusters._layer >= %i && _clusters._layer < %i", layers[i], layers[i+1] ) );

    Utils::TreeToHisto( tree, hname, var2d, momentum_cut && detector_cut && beta_cut && residual_cut, false );
    h->SetTitle(TString());
    h->GetXaxis()->SetTitle( "tan(#beta)" );
    h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );

    cv->cd( i+1 );
    h->Draw();

    auto fname = Form( "f%i", i );
    auto f = new TF1( fname, ::function, -max_beta, max_beta, 2 );
    f->SetParameter( 0, 0 );

    const float r = 0.5*(radius[layers[i]] + radius[layers[i+1]-1]);

    // define the r transformation (cm)
    static constexpr float rmin_tpc = 20;
    static constexpr float rmax_tpc = 78;
    static constexpr float rlength_tpc = rmax_tpc - rmin_tpc;
    static constexpr float deltar_max = 1.5;
    // const float dr = deltar_max*std::sin( M_PI*(r-rmin_tpc)/rlength_tpc );
    const float dr = deltar_max*std::cos( 2*M_PI*(r-rmin_tpc)/rlength_tpc );

    std::cout << "DeltaZ_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta r: " << dr << std::endl;
    f->SetParameter( 1, dr );
    f->SetLineColor(2);
    f->Draw("same");

  }

  pdfDocument.Add( cv );

  return pdfFile;

}
