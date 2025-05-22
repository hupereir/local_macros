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
TString TrackQuality_Tony()
{

  // const TString tag = "_flat_acts_truth_notpc_nodistortion";
  // const TString tag = "_uncorrected-0";
  // const TString tag = "_corrected_module_edge-0";
  // const TString tag = "_corrected-1";
  const TString tag = "Xudong";
  TString postfix;
  TString run_label;


  FileManager fileManager;

  if( false )
  {
    const std::vector<int> runlist = { 53199 };
    // const std::vector<int> runlist = { 53285 };

    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR_CombinedDataReconstruction%s/TrackResiduals-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TrackQuality_Tony - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const std::vector<int> runlist = { 52401 };
    // const std::vector<int> runlist = { 53285 };

    const TString inputFile = "DST/clusters_seeds_52401-0.root_resid.root";
    std::cout << "TrackQuality_Tony - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

//     const TString inputFile = Form( "DST/CONDOR%s/TrackResiduals_flat_acts_truth_notpc_nodistortion_?.root", tag.Data() );
//     std::cout << "TrackQuality_Tony - inputFile: " << inputFile << std::endl;
//     fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TrackQuality_Tony%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "residualtree" );
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
    LogY = 1<<0
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
    { "m_crossing", 100, 0, 100, LogY },
    { "m_pt", 100, 0, 20, LogY },
    { "quality", 100, 0, 500 },
    { "m_ntpc", 50, 0, 50 },
    { "m_nmaps", 5, 0, 5},
    { "m_nintt", 5, 0, 5},
    { "nmms", 5, 0, 5, LogY}
  };

  const TCut trcut(
    "(m_crossing == 1 || m_crossing == 0)"
    "&& m_pt>0.2 "
    "&& quality<100"
    // "&& m_ntpc>20"
    // "&& m_nmaps>2"
    // "&& m_nintt>1"
    // "&& nmms>0"
    );

  for( int i=0; i<vars.size(); ++i )
  {
    const auto& var = vars[i];

    std::cout << "TrackQuality_Tony - var: " << var.m_var << std::endl;

    // canvas
    const auto cvname = Form( "cv_%i", i );
    auto cv = new TCanvas(cvname, "", 900, 900);

    // histogram
    const auto hname = Form( "h_%i", i);
    auto h = new TH1F(hname, "", var.m_nbins, var.m_min, var.m_max);
    Utils::TreeToHisto( tree, hname, var.m_var, trcut, false );
    h->SetTitle("");
    h->GetXaxis()->SetTitle( var.m_var );
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow);
    h->Draw("hist");

    if( var.m_options&LogY) { gPad->SetLogy(true); }

    pdfDocument.Add(cv);
  }

  return pdfFile;

}
