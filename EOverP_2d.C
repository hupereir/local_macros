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

// #include <g4eval_hp/TrackingEvaluator_hp.h>

TString e_over_p_2d()
{

  set_style( false );
  gStyle->SetOptStat(0);

  const TString tag = "_single_electron";
  // const TString tag = "_single_piplus";
  // const TString tag = "_single_proton";

  const TString inputFile = Form( "DST/CONDOR%s/DST_RECO/dst_reco*.root", tag.Data() );
  const TString pdfFile = Form( "Figures/e_over_p_2d%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  std::cout << "e_over_p_2d - entries: " << tree->GetEntries() << std::endl;

  // int entries = tree->GetEntries();
  int entries = 25000;

  // setup branch
  TrackingEvaluator_hp::Container* container = nullptr;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  auto h = new TH2F( "h2d", "h2d", 100, 0, 2, 100, 0, 2 );
  h->GetXaxis()->SetTitle( "E_{EMCAL}/p" );
  h->GetYaxis()->SetTitle( "E_{iHCAL}/E_{EMCAL}" );

  // loop over entries
  for( int ientry = 0; ientry < entries; ++ientry )
  {
    if( !(ientry%1000) ) std::cout << "entry: " << ientry << std::endl;
    tree->GetEntry(ientry);

    for( const auto& track : container->tracks() )
    {

      // check pid
      if( track._pid != 11 ) { continue; }

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

  // make the plot
  auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.15);

  // h->SetMarkerStyle(20);
  h->Draw("SCAT");
  pdfDocument.Add(cv);

  return pdfFile;

}
