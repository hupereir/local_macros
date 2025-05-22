#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>

#include <TPRegexp.h>
#include <TString.h>

#include <fstream>
#include <string>
#include <map>
#include <vector>

R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{
  struct pruner_statistics_t
  {
    uint64_t total_clusters = 0;
    uint64_t removed_clusters = 0;
    double fraction = 0;
  };

  // streamer
  std::ostream& operator << (std::ostream& out, const pruner_statistics_t& statistics )
  {
    out
      << " total: " << statistics.total_clusters
      << " removed: " << statistics.removed_clusters
      << " fraction: " << statistics.fraction;
    return out;
  }
}

//________________________________________________________________________________
pruner_statistics_t ParseFile( const std::string& filename )
{
  std::cout << "ParseFile - filename: " << filename << std::endl;
  TPRegexp regexp_total( "MvtxClusterPruner::End - m_cluster_counter_total: (\\d+)" );
  TPRegexp regexp_removed(
    "MvtxClusterPruner::End - m_cluster_counter_deleted: (\\d+) "
    "fraction: ((?:\\d|.)+)" );

  pruner_statistics_t out;
  std::ifstream in( filename.c_str() );
  std::string line;

  while( std::getline( in, line ) )
  {
    {
      auto matches = regexp_total.MatchS( line.c_str() );
      if( matches->GetLast() >=1 )
      {
        out.total_clusters = static_cast<TObjString*>(matches->At(1))->GetString().Atoi();
      }
    }

    {
      auto matches = regexp_removed.MatchS( line.c_str() );
      if( matches->GetLast() >=1 )
      {
        out.removed_clusters = static_cast<TObjString*>(matches->At(1))->GetString().Atoi();
        out.fraction = static_cast<TObjString*>(matches->At(2))->GetString().Atof();
      }
    }
  }
  return out;
}

//________________________________________________________________________________
TString ParseClusterPrunerStatistics()
{
  const TString tag = "_CombinedDataReconstruction_kfp_test3";
  const TString inputFile = Form( "LOG/kfp_test3/output%s_*_full.txt", tag.Data());
  TString pdfFile( Form( "Figures/ClusterPrunerStatistics%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  std::cout << "ParseTime - loaded " << fileManager.GetNFiles() << " files" << std::endl;

  auto h = new TH1F("h", "", 100, 0, 1 );
  h->GetXaxis()->SetTitle( "Removed clusters fraction" );

  for( const auto& file:fileManager.GetFiles() )
  {
    auto out = ParseFile( file.Data() );
    h->Fill(out.fraction);
  }

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h->SetFillStyle(1001);
  h->SetFillColor(kYellow);
  h->Draw();
  gPad->SetLogy(false);
  pdfDocument.Add(cv);
  return pdfFile;
}
