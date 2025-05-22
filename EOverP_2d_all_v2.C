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

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)


TString e_over_p_2d_all_v2()
{

  set_style( false );
  // gStyle->SetOptStat(0);

  using data_type_pair_t = std::tuple<int, std::string,int>;
  using data_list_t = std::vector<data_type_pair_t>;
  data_list_t data_list =
  {
    { 211, "_single_piplus", kOrange-3 },
    { -211, "_single_piminus", kOrange-4 },
    { 321, "_single_kplus", kCyan-3 },
    { -321, "_single_kminus", kCyan-6 },
    { 2212, "_single_proton", kGreen-2 },
    { -2212, "_single_antiproton", kGreen-3 },
    { 11, "_single_electron", kBlue-4},
    { -11, "_single_positron", kBlue-7 }
  };

  const std::string global_tag = "";

  // pdf document
  const TString pdfFile = Form( "Figures/e_over_p_2d_all_v2%s.pdf", global_tag.c_str());
  PdfDocument pdfDocument( pdfFile );

  // create canvas
  auto cv( new TCanvas( "cv", "cv", 1200, 800 ) );
  TrackingEvaluator_hp::Container* container = nullptr;

  bool first;
  for( const auto& [pid, tag, color]:data_list )
  {

    const TString inputFile = Form( "DST/CONDOR%s/DST_RECO%s/dst_reco*.root", tag.c_str(), global_tag.c_str() );

    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "T" );
    if( !tree ) { continue; }

    // setup branch
    tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

    auto h = new TH2F( Form("h2d_%s",tag.c_str()), "", 50, 0, 2, 50, 0, 2 );
    h->GetXaxis()->SetTitle( "E_{EMCAL}/p" );
    h->GetYaxis()->SetTitle( "E_{iHCAL}/E_{EMCAL}" );

    // int entries = tree->GetEntries();
    int entries = 2500;

    // loop over entries
    for( int ientry = 0; ientry < entries; ++ientry )
    {
      if( !(ientry%1000) ) std::cout << "entry: " << ientry << std::endl;
      tree->GetEntry(ientry);

      for( const auto& track : container->tracks() )
      {

        // check pid
        if( track._pid != pid ) { continue; }

        // momentum
        const auto momentum = track._p;

        double e_emcal = -1;
        double e_ihcal = -1;

        // loop over calo clusters
        for( const auto calo_cluster:track._calo_clusters )
        {
          if( calo_cluster._layer == 1 ) {
            e_emcal = calo_cluster._e;
          } else if( calo_cluster._layer == 2 ) {
            e_ihcal = calo_cluster._e;
          }
        }

        // check energies
        if( e_emcal < 0 || e_ihcal < 0 ) continue;

        // fill histogram
        h->Fill( e_emcal/momentum, e_ihcal/e_emcal );

      }
    }

    h->SetMarkerColor(color);
    h->SetLineColor(color);
    if( first )
    {
      first = false;
      h->Draw("BOX");
    } else {
      h->Draw("BOX SAME");
    }

  }
  pdfDocument.Add(cv);

  return pdfFile;

}
