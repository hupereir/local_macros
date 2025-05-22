#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

namespace
{

  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  template<class T>
    inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  static const std::vector<double> radii = {
    20, 80,
    93.5, 114,
    115, 140.005,
    177.4, 274.317
  };

  void draw_detectors_xy()
  {
    for( const auto& radius:radii )
    {
      auto ellipse = new TEllipse( 0, 0, radius, radius );
      ellipse->SetFillStyle(0);
      ellipse->SetLineColor(1);
      ellipse->Draw();
    }

  }


  void draw_detectors_rz()
  {
    for( const auto& radius:radii )
    {
      auto line = new TLine( -200, radius, 200, radius );
      line->SetLineColor(1);
      line->Draw();
    }

    for( const auto& radius:radii )
    {
      auto line = new TLine( -200, -radius, 200, -radius );
      line->SetLineColor(1);
      line->Draw();
    }
  }

}

//____________________________________________________________________________
void EventDisplay_clusters()
{

  // const TString tag = "_single_electron";
  const TString tag = "_single_piminus";
  const TString inputfile = Form( "DST/CONDOR%s/DST_RECO/dst_reco%s_test_0001.root", tag.Data(), tag.Data() );
  const TString pdffilename = Form( "Figures/EventDisplay_clusters%s.pdf", tag.Data() );

  std::cout << "EventDisplay - inputfile: " << inputfile << std::endl;
  std::cout << "EventDisplay - pdffilename: " << pdffilename << std::endl;

  static const float min_energy_emc = 0.5;
  static const float min_energy_hcalin = 0.0;
  static const float min_energy_hcalout = 0.3;
  // static const float min_energy = 0;

  // pdf output
  PdfDocument pdfDocument( pdffilename );

  // open TFile, load T tree
  FileManager f( inputfile );
  auto tree = f.GetChain( "T" );

  // create container
  TrackingEvaluator_hp::Container* container = nullptr;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  const int entries = 50;
  // const int entries = tree->GetEntries();
  std::cout << "EventDisplay - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    std::cout << "EventDisplay - entry: " << i << std::endl;

    // load entry
    tree->GetEntry(i);
    const auto& clusters = container->clusters();
    const auto& calo_clusters = container->calo_clusters();

    std::cout << "EventDisplay - entry: " << i << " clusters: " << clusters.size() << " calo clusters: " << calo_clusters.size() << std::endl;

    // create canvas
    auto cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    cv->cd(1);
    TH2* h0( new TH2I( "h0", "", 200, -300, 300, 200, -300, 300 ) );
    h0->SetTitle( Form( "event: %i", i ) );
    h0->GetXaxis()->SetTitle( "x (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->SetStats(0);
    h0->Draw();

    draw_detectors_xy();

    std::unique_ptr<TMarker> marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(0.6);
    marker->SetMarkerColor(2);

    // loop over clusters
    for( const auto& cluster:clusters )
    {
      if( cluster._layer != 56 )
      { marker->DrawMarker( cluster._x, cluster._y ); }
    }

    // loop over calo clusters
    for( const auto& cluster:calo_clusters )
    {

      if( cluster._layer == 1 && cluster._e < min_energy_emc ) continue;
      if( cluster._layer == 2 && cluster._e < min_energy_hcalin ) continue;
      if( cluster._layer == 3 && cluster._e < min_energy_hcalout ) continue;
      marker->DrawMarker( cluster._x, cluster._y );
    }

    cv->cd(2);
    TH2* h1( new TH2I( "h1", "", 100, -200, 200, 200, -300, 300 ) );
    h1->GetXaxis()->SetTitle( "z (cm)" );
    h1->GetYaxis()->SetTitle( "r (cm)" );
    h1->SetStats(0);
    h1->Draw();

    draw_detectors_rz();

    for( const auto& cluster:clusters )
    {
      if( cluster._layer != 55 )
      {
        marker->DrawMarker( cluster._z,  (cluster._y > 0) ? cluster._r:-cluster._r );
      }
    }

    for( const auto& cluster:calo_clusters )
    {
      if( cluster._layer == 1 && cluster._e < min_energy_emc ) continue;
      if( cluster._layer == 2 && cluster._e < min_energy_hcalin ) continue;
      if( cluster._layer == 3 && cluster._e < min_energy_hcalout ) continue;
      const auto r = get_r( cluster._x, cluster._y );
      marker->DrawMarker( cluster._z,  (cluster._y > 0) ? r:-r );
    }

    pdfDocument.Add( cv );
  }

}

