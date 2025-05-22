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

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
double fit_function( double* x, double* par )
{
  const double ta = x[0];
  return par[0]+ta*par[1];
}

//____________________________________________________________________________
TString DeltaZBeta_cluster( TString tag = TString() )
{

  set_style( false );

  // mvtx tilts
  std::array<float, 3> tilts = {{ 0.232547, 0.295386, 0.297307 }};
  
  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.01, 0.1, 1.2, 1.2, 1.2, 1.5};
  constexpr float max_beta = 0.9;

  // input files
  if( tag.IsNull() ) tag = "_2k_realistic_full_zhengyun_fix3";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaZBeta_cluster%s.pdf", tag.Data() );

  std::cout << "DeltaZBeta_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZBeta_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
  const TString var3d = Form( "%s:tan( _tracks._clusters._trk_beta ):_tracks._clusters._layer", var.Data() );
  const TCut cluster_cut;
  const TCut momentum_cut;

  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    if( idet != 0 ) continue;

    const TString hname( Form( "deltaz_%i_0", idet ) );
    const TCut layer_cut( Form( "_tracks._clusters._layer==%i", firstLayer[idet]+1 ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut, false );
      max_det_residual[idet] = 5*h1->GetRMS() + std::abs(h1->GetMean() );

    }

  }

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    if( idet != 0 ) continue;

    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH3> h3d( new TH3F( hname, "",
      nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet],
      100, -1, 1,
      100, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var3d, cluster_cut&&momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, (idet == 0) ? 400:800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      h3d->GetXaxis()->SetRange( ilayer+1, ilayer+1 );
      auto h = h3d->Project3D( "zy" );
      h->SetTitle( "" );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "tan(#beta)" );
      h->GetYaxis()->SetTitle( "#Deltaz(cluster-truth) (cm)" );

      cv->cd( ilayer+1 );
      h->Draw();

      if( idet == 0 )
      {
        h->GetYaxis()->SetMaxDigits( 2 );
        h->GetYaxis()->SetTitleOffset( 1.2 );
      }

      // root fit
      auto fname = Form( "f_fit%i", ilayer );
      auto f = new TF1( fname, ::fit_function, -max_beta, max_beta, 2 );
      f->SetParameter( 0, 0 );
      f->SetParameter( 1, 0 );

      h->Fit( f, "0QR" );
      f->SetLineColor( 2 );
      f->Draw("same");

      std::cout << "layer: " << layerIndex << " slope: " << f->GetParameter(1) << " +/- " << f->GetParError(1) << std::endl;
      
      if( layerIndex < 3 ) 
      { std::cout << "layer: " << layerIndex << " slope: " << f->GetParameter(1)/std::cos(tilts[layerIndex]) << " +/- " <<  f->GetParError(1)/std::cos(tilts[layerIndex]) << std::endl; }
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
