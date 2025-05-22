#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

#include <Eigen/Dense>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

void ClusterMultiplicity_micromegas( )
{

//   const TString tag = "_Hijing_Micromegas_50kHz";
//   const TString inputFile = Form( "DST/CONDOR%s/dst_clustereval/3/dst_clustereval*.root", tag.Data() );

  const TString tag = "_Hijing_Micromegas_0_12fm";
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_sHijing_0_12fm/dst_eval_*.root";

  const TString pdfFile = Form( "Figures/ClusterMultiplicity_micromegas%s.pdf", tag.Data() );

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // get number of events
  const auto entries = tree->GetEntries();
  std::cout << "Multiplicity_micromegas - entries: " << entries << std::endl;

  {
    const TString var =
      "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._truth_phi:"
      "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._truth_z";
    const TCut cut = "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._layer==55";
    auto h = new TH2F( "radiograph", "radiograph", 100, -105.5, 105.5, 120, -M_PI, M_PI );

    Utils::TreeToHisto( tree, h->GetName(), var, cut, false );
    h->SetTitle("");
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "#phi (rad)" );
    h->SetMinimum(0);

    auto cv = new TCanvas( "cv0", "cv0", 800, 800 );
    cv->SetRightMargin( 0.2 );
    h->Draw( "colz" );
    pdfDocument.Add( cv );

  }

  {
    const TString var = "DST#EVAL#TrackingEvaluator_hp::Container._events[]._nclusters[56]";
    const TCut cut;

    auto h = new TH1F( "tiles", "tiles", 200, 0, 200 );
    Utils::TreeToHisto( tree, h->GetName(), var, cut, false );
    h->SetTitle("");
    h->GetXaxis()->SetTitle( "clusters/event" );

    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    h->Draw( "h" );
    pdfDocument.Add( cv );

    std::cout << "ClusterMultiplicity_micromegas - clusters/event: " << h->GetMean() << std::endl;

  }


}
