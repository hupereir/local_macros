#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <RootUtil/RootFile.h>


#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
float square(float x) { return x*x; }

//____________________________________________________________________________
float get_r( float x, float y ) { return std::sqrt(square(x)+square(y)); }

//____________________________________________________________________________
float get_eta( float theta )
{ return -std::log(std::tan(theta/2)); }

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString make_run_label( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "run %i", runlist[0] );
  return Form( "runs %i-%i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
TString make_run_postfix( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "_%08i", runlist[0] );
  return Form( "_%08i-%08i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//_____________________________________________
void TrackStatistics_hp()
{

  // const TString tag = "_CombinedDataReconstruction_zf_test";
  const TString tag = "_CombinedDataReconstruction";
  TString postfix;
  TString run_label;

  FileManager fileManager;

  if( true )
  {
    // zero field runs
    const std::vector<int> runlist = { 53285 };
    // const std::vector<int> runlist = { 52077, 52078 };
    for( const auto& runnumber:runlist )
    {
      // const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*.root", tag.Data(), runnumber );
      std::cout << "TrackStatistics_hp - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/TrackResiduals_flat_acts_truth_notpc_nodistortion_?.root", tag.Data() );
    std::cout << "TrackQuality_hp - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  auto tree = fileManager.GetChain( "T" );
  std::cout << "TrackStatistics_hp - tag: " << tag << " label: " << run_label << std::endl;
  std::cout << "TrackStatistics_hp - entries: " << tree->GetEntries() << std::endl;

  const TCut base_trk_cut(
    "(_tracks._crossing == 1 || _tracks._crossing == 0)"
    "&& _tracks._nclusters_tpc>20"
    );

  const TCut match_silicon_cut(
    "_tracks._nclusters_mvtx>=2"
    "&& _tracks._nclusters_intt>=1"
    );

  const TCut match_tpot_cut(
    "_tracks._nclusters_micromegas>=1"
    );

  // get total number of tracks and get silicon matched tracks
  // get total number of tracks that satisfy the cut
  const TString var = "_tracks._pt";
  auto h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut, true );
  auto ntracks = h->GetEntries();
  std::cout << "TrackStatistics_hp - tracks total: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_tpot_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_hp - tracks matched tpot: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_hp - tracks matched silicon: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut && match_tpot_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_hp - tracks matched silicon+tpot: " << ntracks << std::endl;

}
