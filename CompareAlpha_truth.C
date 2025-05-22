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

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString CompareAlpha_truth()
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.003, 0.01, 0.2, 0.2, 0.2, 0.5, 0.5};

  // input files
//   const TString tag = "_acts_truth_notpc_nodistortion";
//   const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco_realistic_micromegas_1?.root", tag.Data() );

  const TString tag = "_genfit_truth_notpc_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco_realistic_micromegas_1?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/CompareAlpha_truth%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/CompareAlpha_truth%s.root", tag.Data() );

  std::cout << "CompareAlpha_truth - inputFile: " << inputFile << std::endl;
  std::cout << "CompareAlpha_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var = "_tracks._clusters[]._trk_alpha";
  const TString var2d = Form( "%s:_tracks._clusters[]._truth_alpha", var.Data() );
  const TCut layer_cut = ("_tracks._clusters[]._layer > 6 && _tracks._clusters[]._layer < 55" );
  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>=2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    );
  
  const double maxalpha= 0.4;
  
  const TString hname = "deltaalpha"; 
  TH2* h = new TH2F( hname, hname, 500, -maxalpha, maxalpha, 500, -maxalpha, maxalpha );
  Utils::TreeToHisto( tree, hname, var2d, layer_cut && momentum_cut && pattern_cut, false );
  
  h->SetTitle("");
  h->GetXaxis()->SetTitle( "tan(#alpha)_{truth}" );
  h->GetYaxis()->SetTitle( "tan(#alpha)_{track}" );
  
  
  const TString cvName = "cv";
  auto cv( new TCanvas( cvName, cvName, 800, 800 ) );

  h->Draw( "colz" );
  pdfDocument.Add( cv );
  
  return pdfFile;

}
