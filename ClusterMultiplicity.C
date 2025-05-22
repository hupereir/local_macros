#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString ClusterMultiplicity()
{

  set_style( false );

//   const TString tag = "_Hijing_Micromegas_50kHz-new";
//   const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_sHijing_0_12fm_50kHz_bkg_0_12fm/dst_eval_*.root";

  const TString tag = "_Hijing_Micromegas-new";
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_sHijing_0_12fm/dst_eval_*.root";

  const TString pdfFile = Form( "Figures/ClusterMultiplicity%s.pdf", tag.Data() );

  std::cout << "ClusterMultiplicity - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  if( !tree ) return TString();
  
  // get number of clusters vs layer
  const TString var = "_clusters._layer";
  auto h = new TH1F( "layers", "", 60, 0, 60 );
  Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );
  h->GetXaxis()->SetTitle( "layer" );

  // number of events
  auto norm = tree->GetEntries();
  std::cout << "Number of events: " << norm << std::endl;

  h->Scale( 1./norm );
  h->SetMinimum(0);

  {
    // create canvas
    const TString cvName = "cv";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    h->Draw();
    pdfDocument.Add( cv.get() );
  }

  // print multiplicity per event for micromegas
  // number of tiles
  static constexpr int ntiles = 4;
  const double cmult[2] =
  {
    h->GetBinContent(56)/ntiles,
    h->GetBinContent(57)/ntiles
  };

  std::cout << "ClusterMultiplicity - layer 55: " << cmult[0] << " clusters/event/tile" << std::endl;
  std::cout << "ClusterMultiplicity - layer 56: " << cmult[1] << " clusters/event/tile" << std::endl;

  // get cluster sizes for mictomegas
  const int ncuts = 2;
  const TCut layercut[] =
  {
    "_clusters._layer == 55",
    "_clusters._layer == 56"
  };

  double csize[2] = {0, 0};
  {
    // create canvas
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 800, 400 ) );
    cv->Divide( 2, 1 );

    for( int icut = 0; icut < 2; ++icut )
    {
      const TString hname = Form( "hclus_%i", icut );
      auto h = new TH1F( hname, "", 30, 0, 30 );
      Utils::TreeToHisto( tree, hname, "_clusters[]._size", layercut[icut], false );
      h->SetTitle("");
      h->GetXaxis()->SetTitle( Form( "csize (layer %i)", icut+55 ) );
      cv->cd(icut+1);
      h->Draw();
      Draw::PutText( 0.6, 0.8, Form( "Mean = %.1f #pm %.1f", h->GetMean(), h->GetMeanError() ) );

      gPad->SetLogy( true );

      csize[icut] = h->GetMean();
      std::cout << "ClusterMultiplicity - layer " << icut+55 << ": cluster size: " << csize[icut] << std::endl;
    }
    pdfDocument.Add( cv.get() );
  }

  // get micromegas cluster distribution
  const TString clus_var = "DST#EVAL#TrackingEvaluator_hp::Container._events._nclusters_micromegas";
  auto hdist = new TH1F( "nclusters_micromegas", "", 1000, 0, 1000 );
  Utils::TreeToHisto( tree, hdist->GetName(), clus_var, TCut(), false );
  hdist->SetTitle( "" );
  hdist->GetXaxis()->SetTitle( "N_{clusters, Micromegas}" );

  {
    const TString cvName = "cv";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    hdist->Draw();
    pdfDocument.Add( cv.get() );
  }

  {
    // clusters per event per tile
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 800, 400 ) );
    cv->Divide( 2, 1 );

    for( int i = 0; i < 2; ++i )
    {
      const TString hname = Form( "hclus_%i", i );
      auto h = new TH1F( hname, "", 30, 0, 30 );
      auto var = Form( "DST#EVAL#TrackingEvaluator_hp::Container._events[]._nclusters[%i]/%i", 55+i, ntiles );
      Utils::TreeToHisto( tree, hname, var, TCut(), false );
      h->SetTitle( "" );
      h->GetXaxis()->SetTitle( Form( "N_{clusters, Micromegas} (layer %i)", 55+i ) );

      cv->cd(i+1);
      h->Draw();
      Draw::PutText( 0.6, 0.8, Form( "Mean = %.1f #pm %.1f", h->GetMean(), h->GetMeanError() ) );

      std::cout << "ClusterMultiplicity - layer " << i+55 << " " << h->GetMean() << ": clusters/event/tile" << std::endl;
    }

    pdfDocument.Add( cv.get() );

  }

  {
    // hits per event per tile
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 800, 400 ) );
    cv->Divide( 2, 1 );

    for( int i = 0; i < 2; ++i )
    {
      const TString hname = Form( "hhits_%i", i );
      auto h = new TH1F( hname, "", 256, 0, 256 );
      auto var = Form( "DST#EVAL#TrackingEvaluator_hp::Container._events[]._nhits[%i]/%i", 55+i, ntiles );
      Utils::TreeToHisto( tree, hname, var, TCut(), false );
      h->SetTitle( "" );
      h->GetXaxis()->SetTitle( Form( "N_{hits, Micromegas} (layer %i)", 55+i ) );

      cv->cd(i+1);
      h->Draw();
      Draw::PutText( 0.6, 0.8, Form( "Mean = %.1f #pm %.1f", h->GetMean(), h->GetMeanError() ) );
      Draw::PutText( 0.6, 0.75, Form( "Occ = %.1f %%", 100*h->GetMean()/256 ) );

      std::cout << "ClusterMultiplicity - layer " << i+55 << " " << h->GetMean() << " hits/event/tile - occupancy: " << 100*h->GetMean()/256 << "%" << std::endl;
    }

    pdfDocument.Add( cv.get() );

  }
  return pdfFile;

}
