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
TString DeltaRPhi_2d()
{

  set_style( false );

  // input files
  // const TString tag = "_TrackReconstruction_corrected-new2";
  const TString tag = "_TrackReconstruction-new2";
  const TString inputFile = Form( "DST/CONDOR%s/TRACKS-00053877-*-full.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaRPhi_2d%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

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
  const TString var2d = Form( "%s:_tracks._clusters[]._trk_r", var.Data() );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "1"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    "&&_tracks._clusters[]._trk_r>0"
    "&&_tracks._clusters[]._layer >=7 && _tracks._clusters[]._layer < 55"
    );

  const TCut zcut("fabs(_tracks._clusters[]._trk_z) < 25");
  Utils::max_entries = 10000;

  auto h = new TH2F( "h", "", 16, 20, 80, 100, -2, 2 );
  h->GetXaxis()->SetTitle("r_{trk} [cm]");
  h->GetYaxis()->SetTitle("r_{trk}#Delta#phi_{clus-trk} [cm]");

  // auto h = new TH2F( "h", "", nLayers_tpc, &tpc_radius[0], 100, -2, 2 );
  Utils::TreeToHisto( tree, h->GetName(), var2d, momentum_cut && pattern_cut && zcut, false );

  auto p = new TProfile( "p", "", 16, 20, 80 );
  Utils::TreeToHisto( tree, p->GetName(), var2d, momentum_cut && pattern_cut && zcut, false );

  auto cv( new TCanvas( "cv", "cv", 900, 600 ) );
  cv->SetRightMargin(0.2);

  h->SetStats(0);
  h->Draw("colz");
  p->SetLineColor(1);
  p->SetMarkerColor(1);
  p->SetMarkerStyle(20);
  p->Draw("same");

  gPad->SetLogz(true);
  gPad->Update();
  Draw::HorizontalLine(gPad, 0)->Draw();

  pdfDocument.Add( cv );

  return pdfFile;
}
