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
TString TrackQuality_hp()
{

  const TString tag = "_CombinedDataReconstruction_test3";
  TString postfix;
  TString run_label;

  FileManager fileManager;

  Utils::max_entries = 10000;

  if( true )
  {
    const std::vector<int> runlist = { 53877 };

    for( const auto& runnumber:runlist )
    {
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

  const bool match_silicons = true;
  const TString pdfFile = match_silicons ?
    Form( "Figures/TrackQuality_hp_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TrackQuality_hp%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );
  std::cout << "entries: " << tree->GetEntries() << std::endl;

  /*
   * variables:
   * - crossing
   * - quality
   * - m_ntpc
   * - m_nmaps
   * - m_nintt
   * nmms
   */

  enum option_t
  {
    None = 0,
    LogY = 1<<0,
    Normalize = 1<<1
  };

  class var_data_t
  {
    public:
    var_data_t( const TString& var, int nbins, float min, float max, int options = 0 ):
      m_var(var),
      m_nbins(nbins),
      m_min(min),
      m_max(max),
      m_options(options)
    {}


    TString m_var;
    int m_nbins=100;
    float m_min=0;
    float m_max=0;
    int m_options = 0;
  };

  const std::vector<var_data_t> vars =
  {
    // { "_tracks._ndf", 100, 0, 100 },
    // { "_tracks._chisquare/_tracks._ndf", 100, 0, 300 },
    // { "_tracks._nclusters_micromegas", 5, 0, 5},
//     { "_tracks._crossing", 100, -100, 500 },
//     { "_tracks._pt", 100, 0, 20, LogY },
    { "_tracks._nclusters_tpc", 50, 0, 50 },
    { "_tracks._nclusters_mvtx", 5, 0, 5 },
    { "_tracks._nclusters_intt", 5, 0, 5 }
  };

  const TCut base_trk_cut(
    "1"
//     "(_tracks._crossing == 1 || _tracks._crossing == 0)"
//    "&& _tracks._pt>0.2"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    );

  const TCut match_tpot_cut( "_tracks._nclusters_micromegas>=1" );

  const TCut match_silicon_cut(
    "_tracks._nclusters_mvtx>0"
    "&& _tracks._nclusters_intt>0"
    );

  // get total number of tracks and get silicon matched tracks
  // get total number of tracks that satisfy the cut
  {
    const TString var = "_tracks._pt";
    auto h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut, true );
    auto ntracks = h->GetEntries();
    std::cout << "TrackQuality_hp - tracks: " << ntracks << std::endl;

    h = Utils::TreeToHisto( tree, "h_pt", var, base_trk_cut && match_silicon_cut, true );
    ntracks = h->GetEntries();
    std::cout << "TrackQuality_hp - tracks matched to silicons: " << ntracks << std::endl;

  }

  const TCut trk_cut = match_silicons ?
    base_trk_cut && match_silicon_cut:
    base_trk_cut;

  for( int i=0; i<vars.size(); ++i )
  {

    const auto& var = vars[i];

    std::cout << "TrackQuality_hp - var: " << var.m_var << std::endl;

    TH1* h_all = nullptr;
    TH1* h_tpot = nullptr;

    {
      // canvas
      const auto cvname = Form( "cv_%i", 3*i );
      auto cv = new TCanvas(cvname, "", 900, 900);

      // histogram
      const auto hname = Form( "h_%i", 3*i);
      auto h = new TH1F(hname, "", var.m_nbins, var.m_min, var.m_max);
      Utils::TreeToHisto( tree, hname, var.m_var, trk_cut, false );
      h->SetTitle("");
      h->GetXaxis()->SetTitle( var.m_var );
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow);

      if( var.m_options&Normalize )
      {
        h->Scale( 1.0/h->GetEntries() );
        h->SetMaximum(1.1);
      } else {
        h->SetMaximum(h->GetMaximum()*1.2);
      }

      h->Draw("hist");

      if( var.m_options&LogY) { gPad->SetLogy(true); }

      Draw::PutText( 0.15, 0.8, Form( "%s", run_label.Data() ));
      pdfDocument.Add(cv);

      h_all = h;

    }

    if( false )
    {
      {
        // canvas
        const auto cvname = Form( "cv_%i", 3*i+1 );
        auto cv = new TCanvas(cvname, "", 900, 900);

        // histogram
        const auto hname = Form( "h_%i", 3*i+1);
        auto h = new TH1F(hname, "", var.m_nbins, var.m_min, var.m_max);
        Utils::TreeToHisto( tree, hname, var.m_var, trk_cut && match_tpot_cut, false );
        h->SetTitle("");
        h->GetXaxis()->SetTitle( var.m_var );
        h->SetFillStyle(1001);
        h->SetFillColor(kYellow);

        if( var.m_options&Normalize )
        {
          h->Scale( 1.0/h->GetEntries() );
          h->SetMaximum(1.1);
        } else {
          h->SetMaximum(h->GetMaximum()*1.2);
        }

        h->Draw("hist");

        if( var.m_options&LogY) { gPad->SetLogy(true); }

        Draw::PutText( 0.15, 0.8, Form( "%s", run_label.Data() ));
        pdfDocument.Add(cv);

        h_tpot = h;

      }

      {
        // canvas
        const auto cvname = Form( "cv_%i", 3*i+2 );
        auto cv = new TCanvas(cvname, "", 900, 900);

        // histogram
        const auto hname = Form( "h_%i", 3*i+2);
        auto h = new TH1F(hname, "", var.m_nbins, var.m_min, var.m_max);
        Utils::DivideHistograms(h_tpot, h_all, h, Utils::ErrorMode::STD );
        h->SetTitle("");
        h->GetXaxis()->SetTitle( var.m_var );
        h->SetFillStyle(1001);
        h->SetFillColor(kYellow);
        h->Draw();

        // if( var.m_options&LogY) { gPad->SetLogy(true); }

        Draw::PutText( 0.15, 0.8, Form( "%s", run_label.Data() ));
        pdfDocument.Add(cv);
      }
    }

  }

  return pdfFile;

}
