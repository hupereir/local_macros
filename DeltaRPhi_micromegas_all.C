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

#include <trackbase/TrkrClusterContainer.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

//____________________________________________________________________________
TString make_run_label( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "run %i", runlist[0] );
  return Form( "runs %i-%i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
TString make_run_postfix( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "_%08i", runlist[0] );
  return Form( "_%08i-%08i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhi_micromegas_all()
{

  set_style( false );

  // initial guess for max residuals
  float max_residual = 1;

  // input files
  const TString tag = "_CombinedDataReconstruction_notpc";

  const TString pdfFile = "Figures/DeltaRPhi_micromegas_all.pdf";
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;
  PdfDocument pdfDocument( pdfFile );

  // variable names
  const TString var( "_tracks._clusters._trk_r*delta_phi(_tracks._clusters._trk_phi - _tracks._clusters._phi)" );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._pt>0.5"
    "&&_tracks._nclusters_tpc>20"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._clusters._trk_r > 0"
    "&& _tracks._clusters._layer == 55"
    "&& _tracks._charge<0"
    );

  auto cv = new TCanvas("cv", "cv", 1200, 800 );
  cv->Divide( 3, 2 );

  int index = 0;
  // for( const int& runnumber : { 53534 })
  for( const int& runnumber : { 53534, 53744, 53756, 53877, 53876, 53630 } )
  {

    FileManager fileManager;

    const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
    std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

    auto run_label = make_run_label( {runnumber} );

    // file manager
    auto tree = fileManager.GetChain( "T" );
    if( !tree ) {
      std::cout << "DeltaRPhi - invalid tree" << std::endl;
      return pdfFile;
    }

    // entries
    std::cout << "DeltaRPhi_raw - entries: " << tree->GetEntries() << std::endl;
    Utils::max_entries = 10000;

    const TString hname( "deltarphi" );
    auto h = new TH1F( hname, "", 100, -max_residual, max_residual );
    Utils::TreeToHisto( tree, hname, var, momentum_cut&&pattern_cut, false );
    h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

    cv->cd( ++index );

    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);

    h->Draw();
    gPad->Update();

    auto line = Draw::VerticalLine(gPad, 0);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();

    Draw::PutText(0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h->GetEntries()));

  }

  pdfDocument.Add( cv );


  return pdfFile;

}
