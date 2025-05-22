#include <RootUtil/FileManager.h>
#include <TChain.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void CheckTrees( const char* inputfiles, const TString& treename = "T", int entries = 0, bool remove_invalid = false )
{
  // file manager
  FileManager fileManager( inputfiles );
  fileManager.CheckTree( treename, entries, remove_invalid );
  return;
}
