#include <RootUtil/FileManager.h>
#include <TChain.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void CheckHistograms( TString inputfiles, int entries = 0, bool remove_invalid = false )
{
  // file manager
  FileManager fileManager( inputfiles );
  fileManager.CheckFiles( remove_invalid );
  return;
}
