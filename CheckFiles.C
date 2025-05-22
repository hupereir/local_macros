#include <RootUtil/FileManager.h>
#include <TChain.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void CheckFiles( const char* inputfiles, bool remove_invalid = false )
{
  // file manager
  FileManager fileManager( inputfiles );
  fileManager.CheckFiles( remove_invalid );
  return;
}
