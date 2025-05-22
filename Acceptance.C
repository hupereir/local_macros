#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>
R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval.so)

#include "LayerDefines.h"

static constexpr bool use_micromegas = false;

//_____________________________________________________________
float get_phi( float phi )
{
  if( phi < 0 ) return phi + 2.*M_PI;
  else if( phi > 2.*M_PI ) return phi - 2.*M_PI;
  else return phi;
} 

//_____________________________________________________________
TString Acceptance()
{
 
  set_style( false );
  
  // input files
  const TString tag = "_Hijing_Micromegas_50kHz_truth_notpc" ;
  // const TString inputFile = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/2/dst_recoeval_sHijing_0-12fm_merged_*.root";
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/2/dst_recoeval_sHijing_0-12fm_merged_001*.root";
  // const TString inputFile = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/2/dst_recoeval_sHijing_0-12fm_merged_000000_000100_truth_notpc.root";

  TString pdfFile( 
    use_micromegas ? 
    Form( "Figures/Acceptance%s_mm_300MeV.pdf", tag.Data() ):
    Form( "Figures/Acceptance%s_300MeV.pdf", tag.Data() )
    );
  PdfDocument pdfDocument( pdfFile );

  // get tree
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return pdfFile;

  if( true )
  {
    auto h2 = new TH2F( "hpr", "hpr", 96, 0, 2*M_PI, nLayers_tpc, &tpc_radius[0] );
    h2->GetXaxis()->SetTitle( "#phi (rad)" );
    h2->GetYaxis()->SetTitle( "r (cm)" );
    
    const TString var = "_tracks._clusters._truth_r:get_phi(_tracks._clusters._truth_phi)";
    TCut cut = 
      // "_tracks._pt > 0.0"
      "_tracks._pt > 0.3"
      "&& _tracks._nclusters_mvtx >= 2"
      "&& _tracks._nclusters_intt >= 2";
    
    if( use_micromegas ) cut = cut && TCut( "_tracks._nclusters_micromegas >= 2" );

    Utils::TreeToHisto( tree, h2->GetName(), var, cut, false );
    h2->SetTitle("");
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->SetRightMargin( 0.24 );
    h2->Draw( "colz" );
    
    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const double phi = (i+1)*M_PI/6;
      Draw::VerticalLine( cv, phi )->Draw();
    }

    pdfDocument.Add( cv );
  }

  if( true )
  {
    auto h2 = new TH2F( "hzr", "hzr", 100, -105.5, 105.5, nLayers_tpc, &tpc_radius[0] );
    h2->GetXaxis()->SetTitle( "z (cm)" );
    h2->GetYaxis()->SetTitle( "r (cm)" );
    
    const TString var = "_tracks._clusters._truth_r:_tracks._clusters._truth_z";
    TCut cut = 
      "_tracks._pt > 0.0"
      // "_tracks._pt > 0.5"
      "&& _tracks._nclusters_mvtx >= 2"
      "&& _tracks._nclusters_intt >= 2";
    
    if( use_micromegas ) cut = cut && TCut( "_tracks._nclusters_micromegas >= 2" );

    Utils::TreeToHisto( tree, h2->GetName(), var, cut, false );
    h2->SetTitle("");
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->SetRightMargin( 0.24 );
    h2->Draw( "colz" );
    
    pdfDocument.Add( cv );
  }
  
  return pdfFile;
}
