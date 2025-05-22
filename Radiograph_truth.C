#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"



template< class T>
T normalize_phi( T phi )
{
  while( phi >= 2*M_PI ) phi -= 2*M_PI;
  while( phi < 0 ) phi += 2*M_PI;
  return phi;
}

namespace
{

  double get_sector_phi( int isec ) { return isec*M_PI/6 + M_PI/12; }
}

//____________________________________________________________________________
void Radiograph_truth()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // input files
  const TString tag = "_realistic_micromegas";
  const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc_distortions-coarse/dst_reco*_?.root";

  // pdf output
  const TString pdfFile = Form( "Figures/Radiograph_truth%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  if( true )
  {
    // y vs x
    const TString var( "_tracks._clusters._truth_y:_tracks._clusters._truth_x" );
    const TString hName = "radiograph";
    auto h2d = new TH2F( hName, "", 1000, -90, 90, 1000, -90, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "x (cm)" );
    h2d->GetYaxis()->SetTitle( "y (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "col" );

    pdfDocument.Add( cv );
  }

  if( false )
  {
    // r vs z
    const TString var( "_tracks._clusters._truth_r:_tracks._clusters._truth_z" );
    const TString hName = "radiograph2";
    auto h2d = new TH2F( hName, "", 500, -110, 110, 500, 0, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "z (cm)" );
    h2d->GetYaxis()->SetTitle( "r (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "col" );

    pdfDocument.Add( cv );
  }

  if( true )
  {
    // r vs phi
    const TString var( "_tracks._clusters._truth_r:normalize_phi(_tracks._clusters._truth_phi)" );
    const TString hName = "radiograph3";
    auto h2d = new TH2F( hName, "", 100, 0, 2.*M_PI, nLayers_tpc, &tpc_radius[0] );
    // auto h2d = new TH2F( hName, "", 100, -M_PI, M_PI, 100, 0, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "phi (rad)" );
    h2d->GetYaxis()->SetTitle( "r (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "colz" );

    gPad->SetRightMargin( 0.18 );

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    { Draw::VerticalLine( cv, get_sector_phi(i)+M_PI/12 )->Draw(); }

    pdfDocument.Add( cv );
  }

  if( true )
  {
    // phi
    const TString var( "normalize_phi(_tracks._clusters[]._phi)" );
    const TCut cut( "_tracks._clusters[]._layer >= 7 && _tracks._clusters[]._layer < 55" );
    const TString hName = "r";
    auto h = new TH1F( hName, "", 500, 0, 2.*M_PI );
    Utils::TreeToHisto( tree, hName, var, cut, kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h->GetXaxis()->SetTitle( "phi (rad)" );
    h->SetTitle( "" );
    h->SetMinimum(0);
    h->Draw();

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    { Draw::VerticalLine( cv, get_sector_phi(i)+M_PI/12 )->Draw(); }

    pdfDocument.Add( cv );
  }

}
