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
TString DeltaRPhi_z_2d()
{

  set_style( false );

  // input files
  // const TString tag = "_TrackReconstruction";
  const TString tag = "_TrackReconstruction_genfit";
  const TString inputFile = Form( "DST/CONDOR%s/TRACKS-00053877-*-full.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_z_2d%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhi_z_2d%s.root", tag.Data() );

  std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = false;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  // variable names
  // const TString var( "_tracks._clusters[]._trk_r*delta_phi(_tracks._clusters[]._trk_phi - _tracks._clusters[]._phi)" );
  const TString var( "_tracks._clusters[]._trk_r*delta_phi(_tracks._clusters[]._phi - _tracks._clusters[]._trk_phi)" );
  const TString var2d = Form( "%s:_tracks._clusters[]._trk_z", var.Data() );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "1"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    "&&_tracks._clusters[]._trk_r>0"
    "&&_tracks._clusters[]._layer >=7 && _tracks._clusters[]._layer < 55"
    );

  Utils::max_entries = 100000;

  auto h = new TH2F( "h", "", 1100, -100, 100, 100, -2, 2 );
  h->GetXaxis()->SetTitle("z_{trk} [cm]");
  h->GetYaxis()->SetTitle("r_{trk}#Delta#phi_{clus-trk}");

  Utils::TreeToHisto( tree, h->GetName(), var2d, momentum_cut && pattern_cut, false );

  auto cv( new TCanvas( "cv", "cv", 900, 600 ) );
  cv->SetRightMargin(0.2);

  h->SetStats(0);
  h->Draw("colz");
  pdfDocument.Add( cv );

  return pdfFile;
}
