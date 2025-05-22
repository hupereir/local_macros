#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
void Radiograph_radial()
{

  set_style( false );

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   const TString tag = "_Hijing_Micromegas_50kHz" ;
//   const TString inputFile = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/dst_recoeval_sHijing_0-12fm_merged_001*.root";

  const TString tag = "_realistic_truth_micromegas_distortions" ;
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );

//   const TString tag = "_realistic_truth_micromegas_nominal" ;
//   const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );

//   const TString tag = "_realistic_truth_micromegas_notpc" ;
//   const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );

  const TString pdfFile = Form( "Figures/Radiograph_radial%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  RootFile rootFile( Form( "Rootfiles/Radiograph_radial%s.root", tag.Data() ) );
      
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  std::cout << "Radiograph_radial - entries: " << tree->GetEntries() << std::endl;
  
  // variable
  const TString var( "_tracks._clusters[]._phi:_tracks._clusters[]._r:_tracks._clusters[]._z" ); 
  
  // track selection cut
  const TCut cut( 
    "_tracks._truth_pt>0.5"
    "&&_tracks._pt>0"
    "&&_tracks._nclusters_mvtx >= 2"
    "&&_tracks._nclusters_intt>=1"
    "&&_tracks._clusters[]._layer[]>=7"
    "&&_tracks._clusters[]._layer[]<55"
    );

  static constexpr int ncut = 2;
  std::array<TCut, ncut> extracut = { TCut(), TCut( "_tracks._nclusters_micromegas>=2" ) };

  for( int i = 0; i < ncut; ++i )
  {
    std::cout << "Radiograph_radial - icut: " << i << std::endl;
    auto hname = Form( "radiograph_%i", i );
    
    // high resolution map
    // auto h = new TH3F( hname, hname, 50, -105, +105, 48, 30, 80, 360, -M_PI, M_PI );
    
    // match axis from ross maps
    auto h = new TH3F( hname, hname, 80, -105.5, +105.5, 16, 20, 78, 36, -M_PI, M_PI );

    // h->GetYaxis()->Set( 48, &tpc_radius[0] );
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "#phi (rad)" );
    
    Utils::TreeToHisto( tree, h->GetName(), var, cut&&extracut[i], false );
    rootFile.Add( h );
  }
  
}
