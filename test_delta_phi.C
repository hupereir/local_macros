#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>

#include <TTree.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

//_________________________________________________________________________
float DeltaPhi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi <= -M_PI ) return phi + 2*M_PI;
  else return phi;
}

void test_delta_phi()
{

  const TString inputFile( "DST/CONDOR_dr_realistic_truth_notpc_noouter/dst_eval_dr_realistic_truth_notpc_noouter_0.root" );
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  tree->Draw("DeltaPhi(_clusters._trk_phi-_clusters._phi)");


}
