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

//_______________________________________________________________
TLine* draw_diagonal( double min, double max )
{
  auto line = new TLine( min, min, max, max);
  line->SetLineStyle( 2 );
  line->SetLineWidth( 1 );
  line->SetLineColor( 1 );
  return line;
}

//_______________________________________________________________
void TpotMatchingCorrelation_test()
{

  set_style( false );

  MicromegasMapping mapping;

//   const TString tag = "_realistic_acts_truth_notpot";
//   TString postfix;
//   const TString inputFile = Form( "DST/CONDOR%s/dst_reco*.root", tag.Data());

  const TString tag = "_CombinedDataReconstruction";
  // const TString tag = "_CombinedDataReconstruction_corrected_module_edge-2";
  TString postfix;
  TString run_label;

  FileManager fileManager;
  if( true )
  {
    // const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };
    // const std::vector<int> runlist = { 51106, 51107, 51109};
    // const std::vector<int> runlist = { 51249 };
    // const std::vector<int> runlist = { 51488};
    // const std::vector<int> runlist = { 52401 };
    // const std::vector<int> runlist = { 53199 };
    const std::vector<int> runlist = { 53285 };
    for( const auto& runnumber:runlist )
    {
      // const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*_corrected.root", tag.Data(), runnumber );
      std::cout << "TpotMatchingCorrelation - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data());
    std::cout << "TpotMatchingCorrelation - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TpotMatchingCorrelation%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  // file manager
  auto tree = fileManager.GetChain( "T" );

  static constexpr int max_x = 15;
  static constexpr int max_y = 30;
  static constexpr int nbins = 100;

  const TCut trk_cut(
    "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._pt>0.2"
    "&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing==1 || DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing==0)"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._ndf > 0"
    "&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._chisquare/DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._ndf)<100"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_tpc>20"
//    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_mvtx>2"
//    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_intt>1"
    );

  gStyle->SetOptStat(111111);
  const bool show_stats = false;

  if(true)
  {
    // local coordinates
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->Divide(4, 4);

    // phi layer
    {
      const TString var = "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_phi._x_local:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._x_local";
      const auto cut = TCut( "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_phi._layer>0 && DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._layer==55" ) && trk_cut;
      const TString var3d = Form( "%s:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_phi._tile",var.Data() );

      const auto hname( "h_phi");
      auto h_3d = new TH3F( hname, "", 8, 0, 8, nbins, -max_x, max_x, nbins, -max_x, max_x );
      Utils::TreeToHisto(tree, h_3d->GetName(), var3d, cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        h_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto h = static_cast<TH2*>(h_3d->Project3D( "zy" ));
        h->SetName( Form( "h_phi_%i", itile ) );
        h->GetXaxis()->SetTitle( "x_{loc} cluster (cm)" );
        h->GetYaxis()->SetTitle( "x_{loc} track (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1);
        h->SetStats(show_stats);
        h->Draw("col");
        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));

        draw_diagonal(-12, 12)->Draw();

      }
    }

    // z layer
    {
      const TString var = "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._y_local:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._y_local";
      const auto cut = TCut( "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._layer>0 && DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._layer==56" ) && trk_cut;
      const TString var3d = Form( "%s:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._tile",var.Data() );

      const auto hname( "h_z");
      auto h_3d = new TH3F( hname, "", 8, 0, 8, nbins, -max_y, max_y, nbins, -max_y, max_y );
      Utils::TreeToHisto(tree, h_3d->GetName(), var3d, cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        h_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto h = static_cast<TH2*>(h_3d->Project3D( "zy" ));
        h->SetName( Form( "h_z_%i", itile ) );

        h->GetXaxis()->SetTitle( "y_{loc} cluster (cm)" );
        h->GetYaxis()->SetTitle( "y_{loc} track (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1+1);
        h->SetStats(show_stats);
        h->Draw("col");
        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));

        draw_diagonal(-25, 25)->Draw();

      }
    }

    pdfDocument.Add(cv);
  }


  if( false )
  {
    // local coordinates
    auto cv = new TCanvas( "cv4", "cv4", 800, 800 );

    // phi layer
    {
      const TString var = "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_phi._x_local:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._x_local";
      const auto cut = TCut( "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_phi._layer>0 && DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._layer==55" ) && trk_cut;

      const auto hname( "h_phi_all");
      auto h = new TH2F( hname, "", nbins, -max_x, max_x, nbins, -max_x, max_x );
      Utils::TreeToHisto(tree, h->GetName(), var, cut, false);

      h->GetXaxis()->SetTitle( "x_{loc} cluster (cm)" );
      h->GetYaxis()->SetTitle( "x_{loc} track (cm)" );
      h->SetTitle("");

      h->SetStats(show_stats);
      h->Draw("col");
      gPad->SetLeftMargin(0.14);
      gPad->SetTopMargin(0.1);
      Draw::PutText( 0.15, 0.92, Form( "#phi views - %s, entries: %.0f", run_label.Data(), h->GetEntries() ));
      draw_diagonal(-12, 12)->Draw();
   }
    pdfDocument.Add(cv);
  }

  if( false )
  {
    // local coordinates
    auto cv = new TCanvas( "cv5", "cv5", 800, 800 );

    // z layer
    {
      const TString var = "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._y_local:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._y_local";
      const auto cut = TCut( "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._layer>0 && DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._clusters[]._layer==56" ) && trk_cut;

      const auto hname( "h_z_all");
      auto h = new TH2F( hname, "", nbins, -max_y, max_y, nbins, -max_y, max_y );
      Utils::TreeToHisto(tree, h->GetName(), var, cut, false);

      h->GetXaxis()->SetTitle( "y_{loc} cluster (cm)" );
      h->GetYaxis()->SetTitle( "y_{loc} track (cm)" );
      h->SetTitle("");

      h->SetStats(show_stats);
      h->Draw("col");
      gPad->SetLeftMargin(0.14);
      gPad->SetTopMargin(0.1);

      Draw::PutText( 0.15, 0.92, Form( "z views - %s, entries: %.0f", run_label.Data(), h->GetEntries() ));
      draw_diagonal(-25, 25)->Draw();
    }
    pdfDocument.Add(cv);
  }

}
