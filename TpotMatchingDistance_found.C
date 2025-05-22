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
void TpotMatchingDistance_found()
{

  set_style( false );

  MicromegasMapping mapping;
  const TString tag = "_CombinedDataReconstruction";
  TString run_label;
  TString postfix;

  FileManager fileManager;
  if( true )
  {
    // const std::vector<int> runlist = { 53534 };
    // const std::vector<int> runlist = { 53744 };
    const std::vector<int> runlist = { 53756 };
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
    Form( "Figures/TpotMatchingDistance_found_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TpotMatchingDistance_found%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );
  std::cout << "TpotMatchingDistance_found - tag: " << tag << " label: " << run_label << std::endl;
  std::cout << "TpotMatchingDistance_found - entries: " << tree->GetEntries() << std::endl;

  static constexpr int max_dx = 10;
  static constexpr int max_dy = 10;

  const TCut base_trk_cut( full_string(
    "1"
//     "&& _tracks._crossing==0"
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

  // get total number of tracks that satisfy the cut
  {
    const TString var = full_string("_tracks._pt");
    auto h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut, true );
    auto ntracks = h->GetEntries();
    std::cout << "TpotMatchingDistance_found - tracks: " << ntracks << std::endl;

    h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut, true );
    ntracks = h->GetEntries();
    std::cout << "TpotMatchingDistance_found - tracks matched: " << ntracks << std::endl;
  }

  {
    // local coordinates
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->Divide(4, 4);

    // phi layer
    {

      const TString var_best = full_string("_tracks._trk_state_phi._x_local-_tracks._best_cluster_phi._x_local");
      const TString var2d_best = Form( full_string("%s:_tracks._trk_state_phi._tile"),var_best.Data() );
      const auto cut_best = TCut( full_string("_tracks._trk_state_phi._layer>0 && _tracks._best_cluster_phi._layer==55")) && trk_cut;
      auto h_2d_best = new TH2F( "h_phi_best", "", 8, 0, 8, 100, -max_dx, max_dx );
      Utils::TreeToHisto(tree, h_2d_best->GetName(), var2d_best, cut_best, false);

      const TString var_found = full_string("_tracks._trk_state_phi._x_local-_tracks._found_cluster_phi._x_local");
      const TString var2d_found = Form( full_string("%s:_tracks._trk_state_phi._tile"),var_found.Data() );
      const auto cut_found = TCut( full_string("_tracks._trk_state_phi._layer>0 && _tracks._found_cluster_phi._layer==55") ) && trk_cut;
      auto h_2d_found = new TH2F( "h_phi_found", "", 8, 0, 8, 100, -max_dx, max_dx );
      Utils::TreeToHisto(tree, h_2d_found->GetName(), var2d_found, cut_found, false);

      std::cout << "TpotMatchingDistance_found - found phi: " << h_2d_found->GetEntries() << std::endl;

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        auto h_best = h_2d_best->ProjectionY( Form( "h_phi_best_%i", itile ), itile+1, itile+1);
        auto h_found = h_2d_found->ProjectionY( Form( "h_phi_found_%i", itile ), itile+1, itile+1);

        for( auto h:{h_best,h_found})
        {
          h->GetXaxis()->SetTitle( "#Deltax_{loc} (track-cluster) (cm)" );
          h->SetTitle("");
        }

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(itile+1);
        h_best->SetFillStyle(1001);
        h_best->SetFillColor(kYellow );
        h_best->Draw();

        h_found->SetFillStyle(1001);
        h_found->SetFillColor(kRed+1 );
        h_found->Draw("same");

        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h_best->GetEntries() ));
      }

    }

    // z layer
    {
      const TString var_best = full_string("_tracks._trk_state_z._y_local-_tracks._best_cluster_z._y_local");
      const TString var2d_best = Form( full_string("%s:_tracks._trk_state_z._tile"),var_best.Data() );
      const auto cut_best = TCut( full_string("_tracks._trk_state_z._layer>0 && _tracks._best_cluster_z._layer == 56" )) && trk_cut;
      auto h_2d_best = new TH2F( "h_z_best", "", 8, 0, 8, 100, -max_dy, max_dy );
      Utils::TreeToHisto(tree, h_2d_best->GetName(), var2d_best, cut_best, false);

      const TString var_found = full_string("_tracks._trk_state_z._y_local-_tracks._found_cluster_z._y_local");
      const TString var2d_found = Form( full_string("%s:_tracks._trk_state_z._tile"),var_found.Data() );
      const auto cut_found = TCut( full_string("_tracks._trk_state_z._layer>0 && _tracks._found_cluster_z._layer == 56")) && trk_cut;
      auto h_2d_found = new TH2F( "h_z_found", "", 8, 0, 8, 100, -max_dy, max_dy );
      Utils::TreeToHisto(tree, h_2d_found->GetName(), var2d_found, cut_found, false);

      std::cout << "TpotMatchingDistance_found - found z: " << h_2d_found->GetEntries() << std::endl;


      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        auto h_best = h_2d_best->ProjectionY( Form( "h_z_best_%i", itile ), itile+1, itile+1);
        auto h_found = h_2d_found->ProjectionY( Form( "h_z_found_%i", itile ), itile+1, itile+1);

        for( auto h:{h_best,h_found})
        {
          h->GetXaxis()->SetTitle( "#Deltay_{loc} (track-cluster) (cm)" );
          h->SetTitle("");
        }

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(itile+MicromegasDefs::m_ntiles+1);
        h_best->SetFillStyle(1001);
        h_best->SetFillColor(kYellow );
        h_best->Draw();

        h_found->SetFillStyle(1001);
        h_found->SetFillColor(kRed+1 );
        h_found->Draw("same");

        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h_best->GetEntries() ));
      }

    }
    pdfDocument.Add(cv);
  }


}
