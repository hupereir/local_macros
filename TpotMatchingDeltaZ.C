#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
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
void TpotMatchingDeltaZ()
{

  set_style( false );

  MicromegasMapping mapping;

  const TString tag = "_CombinedDataReconstruction_corrected_benjamin-new";
  TString postfix;
  TString run_label;

  FileManager fileManager;
  if( true )
  {
    const std::vector<int> runlist = { 53534 };
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
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

  const bool match_silicons = true;

  const TString pdfFile = match_silicons ?
    Form( "Figures/TpotMatchingDeltaZ_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TpotMatchingDeltaZ%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  const TString rootFilename = match_silicons ?
    Form( "Rootfiles/TpotMatchingDeltaZ_matched%s%s.root", tag.Data(),postfix.Data()):
    Form( "Rootfiles/TpotMatchingDeltaZ%s%s.root", tag.Data(),postfix.Data());
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

  {
    // delta z vs z
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);

    static constexpr double max_z = 110;
    static constexpr double max_dz = 10;

    // const TString var = full_string( "_tracks._trk_state_z._z-_tracks._clusters[]._z");
    // const auto cut = TCut( full_string("_tracks._trk_state_z._layer==56 && _tracks._clusters[]._layer==56" )) && trk_cut;
    const TString var = full_string( "_tracks._trk_state_z._z-_tracks._found_cluster_z._z");
    const auto cut = TCut( full_string("_tracks._trk_state_z._layer==56 && _tracks._found_cluster_z._layer==56" )) && trk_cut;

    const TString var2d = Form(full_string("%s:_tracks._trk_state_z._z"), var.Data() );

    TString hname = "h_dz";
    auto h = new TH2F( hname, "", 50, -max_z, max_z, 100, -max_dz, max_dz);
    Utils::TreeToHisto(tree, h->GetName(), var2d, cut, false);
    h->SetTitle("");
    h->GetXaxis()->SetTitle( "z (track) (cm)" );
    h->GetYaxis()->SetTitle( "#Deltaz(track-cluster) (cm)" );
    h->SetStats(0);

    h->Draw("colz" );
    rootFile.Add(h);

    if( true )
    {
      TString pname = "p_dz";
      auto p = h->ProfileX( pname );
      p->SetMarkerStyle(20);
      p->SetMarkerColor(1);
      p->Draw("same");
      rootFile.Add(p);
    }

    // also perform gaussian fit in each x bin
    if( true )
    {
      h->FitSlicesY();
      auto h_mean = static_cast<TH1D*>(gDirectory->Get(Form( "%s_1", h->GetName())));
      auto h_rms = static_cast<TH1D*>(gDirectory->Get(Form( "%s_2", h->GetName())));

      for( int i = 0; i<h_mean->GetNbinsX(); ++i )
      { h_mean->SetBinError(i+1, h_rms->GetBinContent(i+1)/std::sqrt(h->Integral(i+1,i+1,1,h->GetNbinsY()))); }

      h_mean->SetMarkerStyle(20);
      h_mean->SetMarkerColor(2);
      h_mean->SetLineColor(2);
      h_mean->Draw("same");
      rootFile.Add(h_mean);
    }

    gPad->Update();
    Draw::PutText(0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h->GetEntries()));
    Draw::HorizontalLine(gPad, 0)->Draw();

    pdfDocument.Add(cv);
  }
}
