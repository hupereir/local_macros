#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
//
// std::vector<int> runnumbers = {
//     53534,
//     53631,
//     53687,
//     53716,
//     53738,
//     53739,
//     53741,
//     53742,
//     53743,
//     53744,
//     53756,
//     53783,
//     53876,
//     53877,
//     53879
//   };

std::vector<int> runnumbers = {
    53877
  };

void MergeKFPFiles()
{
  const char* tag = "k0s_test1-new";
  for( const auto& runnumber:runnumbers )
  {
    const auto inputFile = Form( "Rootfiles/KFP/%s/KFP-%08i-*-full.root", tag, runnumber );
    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "DecayTree" );
    if(!tree) continue;

    std::cout << "entries: " << tree->GetEntries() << std::endl;

    const auto outputfile = Form( "Rootfiles/KFP/%s/KFP-%08i-merged.root", tag, runnumber );
    tree->Merge(outputfile);
  }
}
