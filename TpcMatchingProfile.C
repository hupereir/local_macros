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
float square( float value ) { return value*value; }

//_______________________________________________________________
float get_r( float x, float y ) { return std::sqrt( square(x)+square(y) ); }

//_______________________________________________________________
void TpcMatchingProfile()
{

  set_style( false );

  MicromegasMapping mapping;

  // const TString tag = "_CombinedDataReconstruction_corrected_notpot";
  const TString tag = "_CombinedDataReconstruction_corrected_notpot_masked";
  TString postfix;
  TString run_label;

  FileManager fileManager;
  if( true )
  {
    const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };
    // const std::vector<int> runlist = { 51106, 51107, 51109};
    // const std::vector<int> runlist = { 51249};
    // const std::vector<int> runlist = { 51488};
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TpcMatchingProfile - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data());
    std::cout << "TpcMatchingProfile - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TpcMatchingProfile%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );
  std::cout << "TpotMatchingEfficiency - entries: " << tree->GetEntries() << std::endl;

  const TCut trk_cut(
    "_tracks._pt>0.2"
    "&& (_tracks._crossing==1 || _tracks._crossing==0)"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    "&& _tracks._nclusters_mvtx>2"
    "&& _tracks._nclusters_intt>1"
    );

  if( true )
  {
    // global coordinates
    auto cv = new TCanvas( "cv", "cv", 1600, 800 );
    cv->Divide(2, 1 );

    // phi layer
    {
      const TString var_rphi = "get_r(_tracks._trk_state_phi._x,_tracks._trk_state_phi._y)*atan2(_tracks._trk_state_phi._y,_tracks._trk_state_phi._x)";
      const TString var_z = "_tracks._trk_state_phi._z";
      const TString var2d = Form( "%s:%s", var_rphi.Data(), var_z.Data() );
      const TCut cut( "_tracks._trk_state_phi._layer==55" );

      TH2* h = new TH2F( "h_phi", "", 100, -120, 120, 100, -220, 20 );
      h->GetXaxis()->SetTitle( "z (cm)" );
      h->GetYaxis()->SetTitle( "r#phi (cm)" );
      Utils::TreeToHisto(tree,h->GetName(), var2d, cut&&trk_cut, false);
      h->SetTitle("");

      cv->cd(1);
      h->Draw("col");
      Draw::PutText( 0.15, 0.85, run_label );

      std::cout << "TpcMatchingProfile - entries: " << h->GetEntries() << std::endl;
    }

    // z layer
    {
      const TString var_rphi = "get_r(_tracks._trk_state_z._x,_tracks._trk_state_z._y)*atan2(_tracks._trk_state_z._y,_tracks._trk_state_z._x)";
      const TString var_z = "_tracks._trk_state_z._z";
      const TString var2d = Form( "%s:%s", var_rphi.Data(), var_z.Data() );
      const TCut cut( "_tracks._trk_state_z._layer==56" );

      auto h = new TH2F( "h_z", "", 100, -150, 150, 100, -250, 50 );
      h->GetXaxis()->SetTitle( "z (cm)" );
      h->GetYaxis()->SetTitle( "r#phi (cm)" );
      Utils::TreeToHisto(tree,h->GetName(), var2d, cut&&trk_cut, false);
      h->SetTitle("");

      cv->cd(2);
      h->Draw("col");
      Draw::PutText( 0.15, 0.85, run_label );

      std::cout << "TpcMatchingProfile - entries: " << h->GetEntries() << std::endl;
    }
    pdfDocument.Add(cv);
  }

  if( true )
  {
    // local coordinates
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->Divide(4, 4);

    // phi layer
    {
      const TString var2d = "_tracks._trk_state_phi._y_local:_tracks._trk_state_phi._x_local";
      const TString var3d = Form( "%s:_tracks._trk_state_phi._tile",var2d.Data() );
      const TCut cut( "_tracks._trk_state_phi._layer>0" );

      const auto hname( "h_phi");
      auto h_3d = new TH3F( hname, "", 8, 0, 8, 100, -30, 30, 100, -30, 30 );
      Utils::TreeToHisto(tree, h_3d->GetName(), var3d, cut&&trk_cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        h_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto h = static_cast<TH2*>(h_3d->Project3D( "zy" ));
        h->SetName( Form( "h_phi_%i", itile ) );
        h->GetXaxis()->SetTitle( "x_{loc} (cm)" );
        h->GetYaxis()->SetTitle( "y_{loc} (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1);
        h->SetStats(0);
        h->Draw("col");
        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));
      }
    }

    // z layer
    {
      const TString var2d = "_tracks._trk_state_z._y_local:_tracks._trk_state_z._x_local";
      const TString var3d = Form( "%s:_tracks._trk_state_z._tile",var2d.Data() );
      const TCut cut( "_tracks._trk_state_z._layer>0" );

      const auto hname( "h_z");
      auto h_3d = new TH3F( hname, "", 8, 0, 8, 100, -30, 30, 100, -30, 30 );
      Utils::TreeToHisto(tree, h_3d->GetName(), var3d, cut&&trk_cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
        h_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto h = static_cast<TH2*>(h_3d->Project3D( "zy" ));
        h->SetName( Form( "h_z_%i", itile ) );
        h->GetXaxis()->SetTitle( "x_{loc} (cm)" );
        h->GetYaxis()->SetTitle( "y_{loc} (cm)" );
        h->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1+1);
        h->SetStats(0);
        h->Draw("col");
        gPad->SetTopMargin(0.1);
        Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));
      }
    }

    pdfDocument.Add(cv);
  }


}
