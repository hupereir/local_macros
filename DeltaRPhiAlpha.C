#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhiAlpha( TString tag = TString() )
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1};

  // input files
  // if( tag.IsNull() ) tag = "_genfit_truth_notpc_nodistortion";
  if( tag.IsNull() ) tag = "_acts_truth_notpc_nodistortion-test2";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco_realistic_micromegas_1?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhiAlpha%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhiAlpha%s.root", tag.Data() );

  std::cout << "DeltaRPhiAlpha - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhiAlpha - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhiAlpha - invalid tree" << std::endl;
    return pdfFile;
  }

  // variable names
  const TString var( "_tracks._clusters._trk_r*delta_phi(_tracks._clusters._trk_phi - _tracks._clusters._phi)" );
  const TString var2d = Form( "%s:tan(_tracks._clusters._trk_alpha)", var.Data() );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>0"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas==2"
    );

  const TCut layer_cut("_tracks._clusters[]._layer>=7 && _tracks._clusters[]._layer<55" );
  
  // histogram
  const TString hname = "h2d";
  const double max_alpha = 0.5;
  const double max_residuals = 0.5;
  auto h2d = new TH2F( hname, hname, 100, -max_alpha, max_alpha, 100, -max_residuals, max_residuals );
  Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut&&layer_cut, false );

  const auto entries = h2d->GetEntries();
  std::cout << "DeltaRPhiAlpha - entries: " << entries << std::endl;
  
  // profile
  const TString pname = "p";
  const TCut residual_cut = Form( "fabs( %s)<%f", var.Data(), max_residuals );
  auto p = new TProfile( pname, pname, 100, -max_alpha, max_alpha );
  Utils::TreeToHisto( tree, pname, var2d, momentum_cut&&pattern_cut&&layer_cut&&residual_cut, false );

  const TString cvName = "cv";
  auto cv( new TCanvas( cvName, cvName, 800, 800 ) );

  h2d->SetTitle("");
  h2d->GetXaxis()->SetTitle( "tan(#alpha)" );
  h2d->GetYaxis()->SetTitle( "r#Delta#phi (track - cluster) (cm)" );
  h2d->Draw();
  
  p->SetLineColor(2);
  p->SetMarkerColor(2);
  p->SetMarkerStyle(20);
  p->Draw("same");
  
  auto line = Draw::HorizontalLine( cv, 0 );
  line->SetLineColor( 4 );
  line->Draw();
  
  pdfDocument.Add( cv );

  return pdfFile;
}
