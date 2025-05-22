#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <micromegas/MicromegasMapping.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

// detector names
std::vector<std::string> detector_names;

MicromegasMapping mapping;

// save detector names
void save_detector_names()
{

  // get detector names that match tile/layer
  for( int itile = 0; itile < 8; ++itile )
  {
    const int layer = 55;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
    const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
    detector_names.push_back( std::move( name.substr(0,3) ) );
  }
}

//_______________________________________________________________
float square( float value ) { return value*value; }

//_______________________________________________________________
float get_r( float x, float y ) { return std::sqrt( square(x)+square(y) ); }

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
void TpotMatchingDeltaZ_tile_fit()
{

  set_style( false );
  save_detector_names();

  const TString tag = "_CombinedDataReconstruction";

  FileManager fileManager;
  const std::vector<int> runlist = { 53534 };
  for( const auto& runnumber:runlist )
  {
    const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
    std::cout << "TpotMatchingCorrelation - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );
  }

  const auto run_label = make_run_label( runlist );
  const auto postfix = make_run_postfix( runlist );

  const bool match_silicons = true;

  const TString pdfFile = match_silicons ?
    Form( "Figures/TpotMatchingDeltaZ_tile_fit_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TpotMatchingDeltaZ_tile_fit%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  const TString rootFilename = match_silicons ?
    Form( "Rootfiles/TpotMatchingDeltaZ_tile_fit_matched%s%s.root", tag.Data(),postfix.Data()):
    Form( "Rootfiles/TpotMatchingDeltaZ_tile_fit%s%s.root", tag.Data(),postfix.Data());
  RootFile rootFile( rootFilename );


  // file manager
  auto tree = fileManager.GetChain( "T" );
  std::cout << "TpotMatchingDeltaZ - Entries: " << tree->GetEntries() << std::endl;

  const TCut base_trk_cut( full_string(
    "1"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    ));

  const TCut match_silicon_cut(full_string(
    "_tracks._nclusters_mvtx>=2"
    "&& _tracks._nclusters_intt>=1"
    ));

  const TCut trk_cut = match_silicons ?
    base_trk_cut && match_silicon_cut:
    base_trk_cut;

  static constexpr double max_z = 110;
  static constexpr double max_dz = 10;

  const TString var = full_string( "_tracks._trk_state_z._z-_tracks._found_cluster_z._z");
  const TString var2d = Form(full_string("%s:_tracks._trk_state_z._z"), var.Data() );
  const TString var3d = Form(full_string("%s:_tracks._trk_state_z._tile"), var2d.Data() );

  const auto cut = TCut( full_string("_tracks._trk_state_z._layer==56 && _tracks._found_cluster_z._layer==56" )) && trk_cut;

  // 3D histogram
  TString hname = "h_dz";
  auto h = new TH3F( hname, "", 8, 0, 8, 50, -max_z, max_z, 100, -max_dz, max_dz);
  // ils::max_entries = 10000;
  Utils::TreeToHisto(tree, h->GetName(), var3d, cut, false);
  h->SetTitle("");
  h->GetXaxis()->SetTitle( "tile number" );
  h->GetYaxis()->SetTitle( "z (track) (cm)" );
  h->GetZaxis()->SetTitle( "#Deltaz(track-cluster) (cm)" );
  h->SetStats(0);

  rootFile.Add(h);

  // canvas
  auto cv = new TCanvas( "cv", "cv", 1200, 800 );
  cv->Divide( 4, 2);
  for( int itile = 0; itile < 8; ++itile )
  {
    cv->cd(itile+1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);

    // project 3D histogram
    h->GetXaxis()->SetRange(itile+1, itile+1);
    auto h2d = static_cast<TH2*>( h->Project3D( "zy" ) );
    const auto hname = Form( "h_%s", detector_names[itile].c_str() );
    h2d->SetName( hname );
    h2d->Draw("colz" );
    h2d->SetStats(0);
    Draw::PutText(0.15, 0.92, Form( "%s %s, entries: %.0f", detector_names[itile].c_str(), run_label.Data(), h->GetEntries()));
    Draw::HorizontalLine(gPad, 0)->Draw();

    rootFile.Add(h2d);
  }
  pdfDocument.Add(cv);
}
