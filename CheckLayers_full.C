#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

TTree* tree = nullptr;

//_________________________________________
TString CheckLayers_full()
{
  
  const TString tag = "_realistic_truth_genfit";
  const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/CheckLayers_full%s.pdf", tag.Data() );

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  tree = fileManager.GetChain( "T" );
  if( !tree ) { 
    std::cout << "CheckLayers_full - invalid tree" << std::endl;
    return pdfFile;
  }
  
  auto container = new TrackingEvaluator_hp::Container;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );


  // loop over entries
  // const int entries = 5e5;
  const int entries = tree->GetEntries();
  std::cout << "CheckLayers_full - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    tree->GetEntry(i);
    if( !(i%100) )
    { std::cout << "CheckLayers_full - entry: " << i << std::endl;}

    // loop over tracks
    for( const auto& track:container->tracks() )
    {
      const auto nclusters1 = track._nclusters_intt;
      const auto nclusters2 = TrackingEvaluator_hp::get_nclusters_intt( track._mask );
      if( nclusters1 != nclusters2 )
      {
        std::cout << "CheckLayers_full -"
          << " entry: " << i 
          << " mvtx: " << nclusters1 
          << " " << nclusters2
          << std::endl;
      }
    }
  }
  
  return pdfFile;

}
