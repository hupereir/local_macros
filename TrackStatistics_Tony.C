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
void TrackStatistics_Tony()
{

  // const TString tag = "_CombinedDataReconstruction_zf_corrected-new";
  const TString tag = "_CombinedDataReconstruction_zf_test";
  TString postfix;
  TString run_label;

  FileManager fileManager;

  {
    // zero field runs
    const std::vector<int> runlist = { 52077, 52078 };
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/TrackResiduals-%08i-*.root", tag.Data(), runnumber );
      std::cout << "TrackStatistics_Tony - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  }

  auto tree = fileManager.GetChain( "residualtree" );
  std::cout << "TrackStatistics_Tony - tag: " << tag << " label: " << run_label << std::endl;
  std::cout << "TrackStatistics_Tony - entries: " << tree->GetEntries() << std::endl;

  const TCut base_trk_cut(
    "(m_crossing == 1 || m_crossing == 0)"
    "&& m_ntpc>20"
    );

  const TCut match_silicon_cut(
    "m_nmaps>=2"
    "&& m_nintt>=1"
    );

  const TCut match_tpot_cut(
    "nmms>0"
    );

  // get total number of tracks and get silicon matched tracks
  // get total number of tracks that satisfy the cut
  const TString var = "m_pt";
  auto h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut, true );
  auto ntracks = h->GetEntries();
  std::cout << "TrackStatistics_Tony - tracks total: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_tpot_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_Tony - tracks matched tpot: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_Tony - tracks matched silicon: " << ntracks << std::endl;

  h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut && match_tpot_cut, true );
  ntracks = h->GetEntries();
  std::cout << "TrackStatistics_Tony - tracks matched silicon+tpot: " << ntracks << std::endl;

}
