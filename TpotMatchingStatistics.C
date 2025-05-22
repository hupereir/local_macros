#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>


#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasDefs.h>

//_______________________________________________________________
float square( float value ) { return value*value; }

//_______________________________________________________________
float get_r( float x, float y ) { return std::sqrt( square(x)+square(y) ); }

//_______________________________________________________________
float get_phi( float x, float y )
{
  return std::atan2(y,x);
}

//_______________________________________________________________
float get_eta( float x, float y, float z )
{ return std::atanh( z/std::sqrt( square(x)+square(y)+square(z) ) ); }

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

//____________________________________________________________________________
TString full_string( const TString& in )
{
  TString out( in );
  return out.ReplaceAll( "_tracks.", "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks." );
}

//_______________________________________________________________
void TpotMatchingStatistics()
{

  set_style( false );

  MicromegasMapping mapping;

  // const TString tag = "_CombinedDataReconstruction_zf_corrected_ppmode_convert_new";
  const TString tag = "_CombinedDataReconstruction_zf_corrected_ppmode_convert_silicon";
  TString run_label;
  TString postfix;

  FileManager fileManager;
  if( true )
  {
    const std::vector<int> runlist = { 52077, 52078 };
    // const std::vector<int> runlist = { 52401 };
    // const std::vector<int> runlist = { 53199 };
    // const std::vector<int> runlist = { 53285 };
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TpotMatchingDistance - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data());
    std::cout << "TpotMatchingDistance - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const bool match_silicons = true;

  const TString pdfFile = match_silicons ?
    Form( "Figures/TpotMatchingStatistics_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TpotMatchingStatistics%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );
  std::cout << "TpotMatchingStatistics - tag: " << tag << " label: " << run_label << std::endl;
  std::cout << "TpotMatchingStatistics - entries: " << tree->GetEntries() << std::endl;

  static constexpr int max_dx = 10;
  static constexpr int max_dy = 10;

  const TCut base_trk_cut( full_string(
    "(_tracks._crossing==1 || _tracks._crossing==0)"
//     "&& _tracks._pt>0.2"
//     "&& _tracks._ndf > 0"
//     "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    ));

  const TCut match_silicon_cut(full_string(
    "_tracks._nclusters_mvtx>=2"
    "&& _tracks._nclusters_intt>=1"
    ));

  const TCut trk_cut = match_silicons ?
    base_trk_cut && match_silicon_cut:
    base_trk_cut;

  // Utils::max_entries = 1000;

  // get total number of tracks that satisfy the cut
  if( false )
  {
    const TString var = full_string("_tracks._pt");
    auto h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut, true );
    auto ntracks = h->GetEntries();
    std::cout << "TpotMatchingStatistics - tracks: " << ntracks << std::endl;

    h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut, true );
    ntracks = h->GetEntries();
    std::cout << "TpotMatchingStatistics - tracks matched: " << ntracks << std::endl;
  }

  // get total number of tracks per tile that satisfy the cut
  if( true )
  {
    const TString var = "_tracks._pt";
    const TString var2d = Form( full_string("%s:_tracks._trk_state_phi._tile"),var.Data() );
    auto h_2d = new TH2F( "h_pt_phi", "", 8, 0, 8, 100, 0, 100 );
    const auto cut = TCut( full_string("_tracks._trk_state_phi._layer>0")) && trk_cut;
    Utils::TreeToHisto(tree, h_2d->GetName(), var2d, cut, false);

    // loop over tiles
    for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
    {
      auto entries = h_2d->ProjectionY( Form( "%s%i", h_2d->GetName(), itile ), itile+1, itile+1)->GetEntries();

      // get detector name
      const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
      const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

      std::cout << "TpotMatchingStatistics - tile: " << itile << " detector: " << detname << " entries: " << entries << std::endl;
    }
    std::cout << std::endl;
  }

  // get total number of tracks per tile that satisfy the cut
  if( false )
  {
    const TString var = "_tracks._pt";
    const TString var2d = Form( full_string("%s:_tracks._trk_state_z._tile"),var.Data() );
    auto h_2d = new TH2F( "h_pt_z", "", 8, 0, 8, 100, 0, 100 );
    const auto cut = TCut( full_string("_tracks._trk_state_z._layer>0")) && trk_cut;
    Utils::TreeToHisto(tree, h_2d->GetName(), var2d, cut, false);

    // loop over tiles
    for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
    {
      auto entries = h_2d->ProjectionY( Form( "%s%i", h_2d->GetName(), itile ), itile+1, itile+1)->GetEntries();

      // get detector name
      const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
      const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

      std::cout << "TpotMatchingStatistics - tile: " << itile << " detector: " << detname << " entries: " << entries << std::endl;
    }
    std::cout << std::endl;
  }

  // extrapolated tracks to TPOT
  if( false )
  {
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    const TString var = full_string("get_phi(_tracks._trk_state_phi._x,_tracks._trk_state_phi._y):_tracks._trk_state_phi._z");
    auto h_2d = new TH2F( "h_track", "", 100, -110, 110, 100, -3.2, 0 );
    const auto cut = TCut( full_string("_tracks._trk_state_phi._layer>0")) && trk_cut;
    Utils::TreeToHisto(tree, h_2d->GetName(), var, cut, false);

    h_2d->GetXaxis()->SetTitle( "z (cm)");
    h_2d->GetYaxis()->SetTitle( "#phi (rad)");
    h_2d->Draw("colz");
    h_2d->SetStats(0);
    cv->SetRightMargin(0.15);
    cv->SetLeftMargin(0.15);
    h_2d->GetYaxis()->SetTitleOffset(1.5);

    gPad->SetTopMargin(0.1);
    Draw::PutText( 0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h_2d->GetEntries() ));

    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    const TString var = full_string("get_phi(_tracks._trk_state_phi._x,_tracks._trk_state_phi._y):get_eta(_tracks._trk_state_phi._x,_tracks._trk_state_phi._y,_tracks._trk_state_phi._z)");
    auto h_2d = new TH2F( "h_track", "", 100, -1.2, 1.2, 100, -3.2, 0 );
    const auto cut = TCut( full_string("_tracks._trk_state_phi._layer>0")) && trk_cut;
    Utils::TreeToHisto(tree, h_2d->GetName(), var, cut, false);

    h_2d->GetXaxis()->SetTitle( "#eta");
    h_2d->GetYaxis()->SetTitle( "#phi (rad)");
    h_2d->Draw("colz");
    h_2d->SetStats(0);
    cv->SetRightMargin(0.15);
    cv->SetLeftMargin(0.15);
    h_2d->GetYaxis()->SetTitleOffset(1.5);

    gPad->SetTopMargin(0.1);
    Draw::PutText( 0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h_2d->GetEntries() ));

    pdfDocument.Add(cv);
  }

}
