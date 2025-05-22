#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;

//____________________________________________________________________________
TString DeltaZ_cluster_2d( TString tag = TString() )
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.003, 5, 1.2, 1.2, 1.2, 0.1, 0.1};

  if( tag.IsNull() ) tag = "_truth_old";
  // if( tag.IsNull() ) tag = "_truth_orig";
  // if( tag.IsNull() ) tag = "_truth_new4";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/DeltaZ_cluster_2d%s.pdf", tag.Data() );

  std::cout << "DeltaZ_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
  const TString var3d = Form( "%s:_tracks._clusters._layer:_tracks._clusters._truth_z", var.Data() );
  const TCut momentum_cut;
  const TCut cluster_cut;
  
  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // skip detector 5, which is phi segmented micromegas
    if( idet == 5 ) continue;

    const TString hname( Form( "deltarphi_%i_0", idet ) );

    const TCut layer_cut( Form( "_tracks._clusters._layer==%i", firstLayer[idet] ) );

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

    std::cout << "DeltaZ_cluster - detector " << idet << " max_det_residual: " << max_det_residual[idet] << std::endl;
    
    // skip detector 5, which is phi segmented micromegas
    if( idet == 5 ) continue;

    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH3> h3d( new TH3F( hname, "", 
      100, m_zmin, m_zmax,
      nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 
      100, -max_det_residual[idet], max_det_residual[idet] ) );
    Utils::TreeToHisto( tree, hname, var3d, cluster_cut&&momentum_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, (idet == 0 ) ? 400:800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      h3d->GetYaxis()->SetRange( ilayer+1, ilayer+1 );
      auto h( h3d->Project3D( "zx" ) );
      h->SetName( hname );
      h->SetTitle( hname );

      h->GetXaxis()->SetTitle( "z_{clus} (cm)" );

      h->GetYaxis()->SetTitle( "#Deltaz_{cluster-truth} (cm)" );
      h->GetYaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
     
      cv->cd( ilayer+1 );
      h->Draw();

      // also project into profile
//       auto profile = h3d->Project3DProfile( "zx" );
      auto profile = static_cast<TH2*>(h)->ProfileX();
      profile->SetLineColor(2);
      profile->SetMarkerColor(2);
      profile->Draw("same");
      
      // fit
      if( false )
      {
        profile->Fit( "pol1", "0Q" );
        auto f = profile->GetFunction( "pol1" );
        f->SetLineColor(2);
        f->Draw("same" );
      }
      
      if( idet == 0 )
      {
        h->GetYaxis()->SetMaxDigits( 2 );
        gPad->SetRightMargin( 0.12 );
      }

      // draw vertical line at zero
      gPad->Update();
      Draw::HorizontalLine( gPad, 0 )->Draw();

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
