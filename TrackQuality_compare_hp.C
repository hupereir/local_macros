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
TString TrackQuality_compare_hp()
{

  // const TString tag = "_flat_acts_truth_notpc_nodistortion";
  const TString tag = "_CombinedDataReconstruction_corrected-2";
  TString postfix;
  TString run_label;


  FileManager fileManager;

  if( true )
  {
    // const std::vector<int> runlist = { 52401 };
    // const std::vector<int> runlist = { 53199 };
    const std::vector<int> runlist = { 53285 };

    for( const auto& runnumber:runlist )
    {
      // const TString inputFile = Form( "DST/CONDOR_CombinedDataReconstruction%s/dst_eval-%08i-*.root", tag.Data(), runnumber );
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TrackQuality_hp - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/TrackResiduals_flat_acts_truth_notpc_nodistortion_?.root", tag.Data() );
    std::cout << "TrackQuality_hp - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TrackQuality_compare_hp%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );
  std::cout << "entries: " << tree->GetEntries() << std::endl;

  const TCut trcut_1(
    "(_tracks._crossing == 1 || _tracks._crossing == 0)"
    "&& _tracks._pt>0.3"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    "&& _tracks._nclusters_mvtx>=2"
    "&& _tracks._nclusters_intt>=1"
    // "&& _tracks._nclusters_micromegas>0"
    );

  const TString var_1 = "_tracks._nclusters_micromegas";

  const TCut trcut_2(
    "(DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._crossing == 1 || DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._crossing == 0)"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._pt>0.2"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._ndf > 0"
    "&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._chisquare/DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._ndf)<100"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._nclusters_tpc>20"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._nclusters_mvtx>=2"
    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._nclusters_intt>=1"
    // "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._nclusters_micromegas>0"
    );

  const TString var_2 = "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks.._nclusters_micromegas";

  {
    // canvas
    const auto cvname( "cv_0");
    auto cv = new TCanvas(cvname, "", 900, 900);

    // histograms
    {
      const auto hname1("h_1");
      auto h1 = new TH1F(hname1, "",  5, 0, 5 );
      Utils::TreeToHisto( tree, hname1, var_1, trcut_1, false );

      h1->SetTitle("");
      h1->GetXaxis()->SetTitle( var_1 );
      h1->SetFillStyle(3003);
      h1->SetFillColor(kYellow);
      h1->Draw("hist");
      h1->SetMaximum(h1->GetMaximum()*1.2);
    }

    {
      const auto hname2("h_2");
      auto h2 = new TH1F(hname2, "",  5, 0, 5 );
      Utils::TreeToHisto( tree, hname2, var_2, trcut_2, false );

      h2->SetTitle("");
      h2->GetXaxis()->SetTitle( var_2 );
      h2->SetFillStyle(3003);
      h2->SetFillColor(kRed+1);
      h2->Draw("hist same");
    }

    gPad->SetLogy(true);

    pdfDocument.Add(cv);
  }

  return pdfFile;

}
