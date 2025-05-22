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

// radial selection
constexpr double max_residual = 0.5;
constexpr double radius_rec = 43.56;

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;


//____________________________________________________________________________
TString DeltaZ_distribution_cluster( TString tag = TString() )
{

  set_style( false );

//   if( tag.IsNull() ) tag = "_acts_full_newgeom";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  
  if( tag.IsNull() ) tag = "_realistic_micromegas";
  const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc/dst_reco*.root";

  const TString pdfFile = Form( "Figures/DeltaZ_distribution_cluster%s.pdf", tag.Data() );

  std::cout << "DeltaZ_distribution_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_distribution_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // find matching layer
  int ilayer = 0;
  for( ; ilayer < nLayersTotal && radius[ilayer] < radius_rec; ++ilayer ) {}
  --ilayer;

  std::cout << "DeltaZ_distribution_cluster - layer: " << ilayer << std::endl;
  
  // variable names
  // const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
  const TString var( "_tracks._clusters._z - _tracks._clusters._trk_z" );
  const TString var2d = Form( "%s:_tracks._clusters._z", var.Data() );
  const TCut momentum_cut ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>0"
    "&&_tracks._nclusters_intt>=2"
    );

  const TCut layer_cut( Form( "_tracks._clusters._layer == %i", ilayer ) );

  // project
  const TString hname = "hdz_z";
  auto h = new TH2F( hname, "", 200, m_zmin, m_zmax, 100, -max_residual, max_residual );
  Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&layer_cut&&pattern_cut, false );

  // also create profile
  const TString pname = "pdz_z";
  auto p = new TProfile( pname, "", 200, m_zmin, m_zmax );
  Utils::TreeToHisto( tree, pname, var2d, momentum_cut&&layer_cut&&pattern_cut, false );
  
  {
    auto cv( new TCanvas( "cv1", "cv1", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );
    cv->SetRightMargin( 0.16 );
    
    h->Draw( "colz" );
    
    p->SetLineColor(1);
    p->SetMarkerColor(1);
    p->Draw( "same" );
    
    pdfDocument.Add( cv );
  }

 
  {
    auto cv( new TCanvas( "cv2", "cv2", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    p->SetMaximum( 0.06 );
    p->SetMinimum( -0.06 );
    p->Draw();
    
    Draw::HorizontalLine( cv, 0 )->Draw();
    
    pdfDocument.Add( cv );
  }
  
  
  return pdfFile;

}
