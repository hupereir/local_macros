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
TString DeltaRPhi_micromegas_raw()
{

  set_style( false );

  // initial guess for max residuals
  float max_residual = 3;

  // input files
  // const TString tag = "_CombinedDataReconstruction_corrected_benjamin";
  const TString tag = "_CombinedDataReconstruction_corrected_benjamin-new";
  // const TString tag = "_CombinedDataReconstruction";
  TString run_label;
  TString postfix;

  FileManager fileManager;

  {
    const std::vector<int> runlist = { 53534 };
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );
  }

  const TString pdfFile = Form( "Figures/DeltaRPhi_micromegas_raw%s%s.pdf", tag.Data(), postfix.Data() );
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;
  PdfDocument pdfDocument( pdfFile );


  // file manager
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  // entries
  std::cout << "DeltaRPhi_raw - entries: " << tree->GetEntries() << std::endl;
  // Utils::max_entries = 10000;

  // variable names
  const TString var( "_tracks._clusters._trk_r*delta_phi(_tracks._clusters._trk_phi - _tracks._clusters._phi)" );
  const TString var2d = Form( "%s:_tracks._clusters._tileid", var.Data() );

  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "1"
//    "&&_tracks._pt>0.5"
//    "&&_tracks._crossing == 0"
    "&&_tracks._nclusters_tpc>20"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
//     "&& _tracks._ndf > 0"
//     "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._clusters._trk_r > 0"
    "&& _tracks._clusters._layer == 55"
    );

  const TString hname( "deltarphi" );
  auto h2d = new TH2F( hname, "",8,0,8, 100, -max_residual, max_residual );
  Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut, false );
  h2d->GetXaxis()->SetTitle( "tile ID" );
  h2d->GetYaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

  {
    // create canvas
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    h2d->Draw("colz");
    pdfDocument.Add( cv );
  }


  {
    auto cv( new TCanvas( "cv2", "cv2", 1200, 800 ) );
    cv->Divide( 4, 2 );
    for( int itile = 0; itile<8; ++itile )
    {
      const auto hname = Form( "h_%i", itile );
      auto h = h2d->ProjectionY(hname, itile+1, itile+1);

      cv->cd(itile+1);
      h->SetMaximum( 1.2*h->GetMaximum() );
      h->Draw();
    }
    pdfDocument.Add( cv );
  }


  return pdfFile;

}
