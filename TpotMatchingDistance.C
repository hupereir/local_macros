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
  return out.ReplaceAll( "_tracks.", "_tracks." );
}

//_______________________________________________________________
void TpotMatchingDistance()
{

  set_style( false );

  MicromegasMapping mapping;

  const TString tag = "_CombinedDataReconstruction_corrected-2";
  TString run_label;
  TString postfix;

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

  const TString pdfFile = Form( "Figures/TpotMatchingDistance%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );

  static constexpr int max_dx = 10;
  static constexpr int max_dy = 10;

  const TCut trk_cut(full_string(
    "_tracks._pt>0.2"
    "&& (_tracks._crossing==1 || _tracks._crossing==0)"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
//    "&& _tracks._nclusters_mvtx>2"
//    "&& _tracks._nclusters_intt>1"
    ));

  {
    // local coordinates
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->Divide(4, 4);

    // phi layer
    {
      const TString var = full_string("_tracks._trk_state_phi._x_local-_tracks._clusters[]._x_local");
      const TString var2d = Form( full_string("%s:_tracks._trk_state_phi._tile",var.Data() ));
      const auto cut = TCut( full_string("_tracks._trk_state_phi._layer>0 && _tracks._clusters[]._layer==55" )) && trk_cut;

      const auto hname( "h_phi");
      auto h_2d = new TH2F( hname, "", 8, 0, 8, 100, -max_dx, max_dx );
      Utils::TreeToHisto(tree, h_2d->GetName(), var2d, cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        auto h = h_2d->ProjectionY( Form( "h_phi_%i", itile ), itile+1, itile+1);

        h->GetXaxis()->SetTitle( "#Deltax_{loc} (track-cluster) (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1);
        h->SetFillStyle(1001);
        h->SetFillColor(kYellow );
        h->Draw();

        h->Fit( "gaus", "0");
        auto f = h->GetFunction("gaus");
        f->SetLineColor(2);
        f->Draw("same");

        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));
      }

    }

    // z layer
    {
      const TString var = full_string("_tracks._trk_state_z._y_local-_tracks._clusters[]._y_local");
      const TString var2d = Form( full_string("%s:_tracks._trk_state_z._tile",var.Data()) );
      const auto cut = TCut( full_string("_tracks._trk_state_z._layer>0 && _tracks._clusters[]._layer == 56") ) && trk_cut;

      const auto hname( "h_z");
      auto h_2d = new TH2F( hname, "", 8, 0, 8, 100, -max_dy, max_dy );
      Utils::TreeToHisto(tree, h_2d->GetName(), var2d, cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        auto h = h_2d->ProjectionY( Form( "h_z_%i", itile ), itile+1, itile+1);

        h->GetXaxis()->SetTitle( "#Deltay_{loc} (track-cluster) (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1+1);
        h->SetFillStyle(1001);
        h->SetFillColor(kYellow );
        h->Draw();
        h->Fit( "gaus", "0");
        auto f = h->GetFunction("gaus");
        f->SetLineColor(2);
        f->Draw("same");

        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));
      }

    }
    pdfDocument.Add(cv);
  }


}
