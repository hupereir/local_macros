#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;

//__________________________________________________________________________
double get_phi( double phi )
{ return phi >= 0 ? phi: (phi+2.*M_PI); }

//__________________________________________________________________________
TString ScanClusters()
{
  set_style( false );
  // const TString tag = "_flat_truth_micromegas_nominal";
  const TString tag = "_flat_truth_micromegas_corrected_mm-coarse_extrapolated-new2";

  // const TString inputFile = Form( "DST/CONDOR%s/dst_eval*_1?.root", tag.Data() );
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data() );
  const TString pdfFile = Form( "Figures/ScanClusters%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // variable
  const TString var( "_clusters._z:_clusters._layer:get_phi(_clusters._phi)" );

  // load input
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // cluster scan
  auto h = new TH3F( "scan", "scan",
    100, 0, 2.*M_PI,
    nLayers_tpc, firstLayer_tpc, firstLayer_tpc + nLayers_tpc,
    100, -110, 110 );

  Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );
  h->SetTitle( "" );
  h->GetXaxis()->SetTitle( "#phi (rad)" );
  h->GetYaxis()->SetTitle( "layer" );
  h->GetZaxis()->SetTitle( "z (cm)" );

  // project on layers
  for( int ilayer = 0; ilayer < nLayers_tpc; ++ilayer )
  {
    int layerIndex = ilayer+firstLayer_tpc;
    const auto hname = Form( "%s_%i", h->GetName(), layerIndex );

    // layer range
    h->GetYaxis()->SetRange( ilayer+1, ilayer+1 );
    const auto hlayer( h->Project3D( "zx" ) );
    hlayer->SetName( hname );
    hlayer->SetTitle("");

    const auto cvname = Form( "cv_%i", ilayer );
    auto cv = new TCanvas( cvname, cvname, 800, 800 );
    cv->SetLeftMargin( 0.14 );
    cv->SetRightMargin( 0.20 );

    hlayer->Draw( "colz" );
    hlayer->GetYaxis()->SetTitleOffset( 1.35 );

    Draw::PutText( 0.2, 0.8, Form( "layer: %i - r = %.2f, entries=%.0f",
      layerIndex,
      radius[layerIndex],
      hlayer->GetEntries() ) );

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const double phi = (i+1)*M_PI/6;
      Draw::VerticalLine( cv, phi )->Draw();
    }

    {
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }

    pdfDocument.Add( cv );

  }

  return pdfFile;

}
