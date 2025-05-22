#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void CountEvents( TString tag = TString() )
{
  FileManager fileManager( "DST/CONDOR_Hijing_Micromegas/G4Hits_sHijing_0-12fm_merged_*.root" );
  auto tree = fileManager.GetChain( "T" );
  std::cout << "CountEvents - entries: " << tree->GetEntries() << std::endl;
  
}
