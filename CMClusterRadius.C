#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TEllipse.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval.so)

//____________________________________________________________________________
void CMClusterRadius()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   // input files
//   const TString tag = "_single_directlasers_test_no_distortion";
//   const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );

  // input files
  const TString tag = "_centralmembrane-nominal";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  // pdf output
  const TString pdfFile = Form( "Figures/CMClusterRadius%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  TH1* h_radius_g4hits = new TH1F( "h_radius_g4hits", "h_radius_g4hits", 500, 20, 80 );
  Utils::TreeToHisto(tree, "h_radius_g4hits", "_g4hits._r", "_g4hits._z>0", false );
  
  TH1* h_radius_clusters = new TH1F( "h_radius_clusters", "h_radius_clusters", 500, 20, 80 );
  Utils::TreeToHisto(tree, "h_radius_clusters", "_clusters._r", "_clusters._z>0", false );

  TH1* h_radius_cm_clusters = new TH1F( "h_radius_cm_clusters", "h_radius_cm_clusters", 500, 20, 80 );
  Utils::TreeToHisto(tree, "h_radius_cm_clusters", "_cm_clusters._r", "_cm_clusters._z>0", false );
  
  const auto maximum = std::max( {h_radius_g4hits->GetMaximum(), h_radius_clusters->GetMaximum(), h_radius_cm_clusters->GetMaximum() } );

  for( auto& h:{h_radius_g4hits, h_radius_clusters, h_radius_cm_clusters} )
  { 
    h->SetMaximum( 1.2*maximum );
    h->GetXaxis()->SetTitle( "r (cm)" ); 
  }
    
  auto cv = new TCanvas("cv", "cv", 900, 900 );

  h_radius_g4hits->SetLineColor( 1 );
  h_radius_g4hits->SetFillColor( 1 );
  h_radius_g4hits->Draw();
  
  h_radius_clusters->SetLineColor( 2 );
  h_radius_clusters->SetFillColor( 2 );
  h_radius_clusters->Draw( "same" );
  
  h_radius_cm_clusters->SetLineColor( 4 );
  h_radius_cm_clusters->SetFillColor( 4 );
  h_radius_cm_clusters->Draw( "same" );
  
  pdfDocument.Add( cv );
}
