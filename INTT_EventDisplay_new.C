#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <g4eval_hp/MicromegasGeometryContainer.h>
#include <g4eval_hp/TrackingEvaluator_hp.h>

//! lines
using line_t = std::array<double,4>;
std::ostream& operator << (std::ostream& out, const line_t& line )
{
  out << "{ " << line[0] << ", " << line[1] << ", " << line[2] << ", " << line[3] << "}";
  return out;
}

namespace
{

  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  template<class T>
    inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  // geometry container
  std::unique_ptr<MicromegasGeometryContainer> m_geometry_container;

  // load geometry from file
  void load_geometry( const std::string& filename )
  {
    auto inputfile( TFile::Open( filename.c_str() ) );
    m_geometry_container.reset( dynamic_cast<MicromegasGeometryContainer*>( inputfile->Get( "MicromegasGeometryContainer" ) ) );
    m_geometry_container->identify();
  }

  /// intt geometry
  std::vector<line_t> intt_segments;

  //_______________________________________________
  void load_intt_geometry( const TString& filename )
  {

    /// running intt segment
    intt_segments.clear();
    auto intt_segment = std::make_unique<line_t>();

    auto tfile = std::make_unique<TFile>( filename, "READ" );
    auto tree = static_cast<TTree*>( tfile->Get("tree"));

    tree->SetBranchAddress( "x1", &(*intt_segment.get())[0] );
    tree->SetBranchAddress( "y1", &(*intt_segment.get())[1] );
    tree->SetBranchAddress( "x2", &(*intt_segment.get())[2] );
    tree->SetBranchAddress( "y2", &(*intt_segment.get())[3] );

    const auto entries = tree->GetEntries();
    for( int i = 0; i < entries; ++i )
    {
      tree->GetEntry( i );
      intt_segments.push_back( *intt_segment );
    }

    std::cout << "load_intt_geometry - intt_segments size: " << intt_segments.size() << std::endl;

  }

  //____________________________________________________________________________
  void draw_intt_xy()
  {
    auto line = std::make_unique<TLine>();
    line->SetLineColor(1);
    for( const auto& [x1,y1,x2,y2]:intt_segments )
    { line->DrawLine( x1/10, y1/10, x2/10, y2/10 ); }
  }

  //____________________________________________________________________________
  void draw_intt_rz()
  {
    for( const double& r:{-7.7, -10.3, 7.7, 10.3} )
    {
      auto line = new TLine( -25, r, 25, r );
      line->SetLineColor(1);
      line->Draw();
    }
  }

  //_______________________________________________
  TVector3 get_strip_center( int layer, int tile, int strip )
  {
    // update cluster position
    const auto strip_begin = m_geometry_container->get_strip_begin( layer, tile, strip );
    const auto strip_end = m_geometry_container->get_strip_end( layer, tile, strip );
    return 0.5*(strip_begin + strip_end);
  }

  //____________________________________________________________________________
  void draw_tpot_xy()
  {

    // loop over tiles
    for( int tile = 0; tile < 8; ++tile )
    {
      const auto first_strip_center = get_strip_center( 55, tile, 0 );
      const auto last_strip_center = get_strip_center( 55, tile, 255 );
      auto line = new TLine( first_strip_center.x(), first_strip_center.y(), last_strip_center.x(), last_strip_center.y() );
      line->SetLineColor(1);
      line->Draw();
    }
  }

  //____________________________________________________________________________
  void draw_tpot_rz()
  {

    // loop over tiles
    for( int tile = 0; tile < 8; ++tile )
    {
      const auto first_strip_center = get_strip_center( 56, tile, 0 );
      const auto last_strip_center = get_strip_center( 56, tile, 255 );
      auto line = new TLine(
        first_strip_center.z(), -get_r(first_strip_center.x(), first_strip_center.y()),
        last_strip_center.z(), -get_r(last_strip_center.x(), last_strip_center.y()) );
      line->SetLineColor(1);
      line->Draw();
    }
  }

}

//____________________________________________________________________________
void INTT_EventDisplay_new( const int runnumber = 25926 )
{

  const TString inputfile = Form( "DST/CONDOR_CombinedDataReconstruction/dst_eval-%08i-0005.root", runnumber );

  const bool swap_tpc_z = false;
  const TString pdffilename = swap_tpc_z ?
    Form( "Figures/INTT_EventDisplay_new-%08i-0005_swapped.pdf", runnumber ):
    Form( "Figures/INTT_EventDisplay_new-%08i-0005.pdf", runnumber );

  // load intt geometry
  load_intt_geometry( "INTT_TLine_geo.root" );

  // load geometry
  load_geometry( "micromegas_geometry.root" );

  // pdf output
  PdfDocument pdfDocument( pdffilename );

  // open TFile, load T tree
  FileManager f( inputfile );
  auto tree = f.GetChain( "T" );

  // create container
  TrackingEvaluator_hp::Container* container = nullptr;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  // const int entries = 5e5;
  const int entries = tree->GetEntries();
  std::cout << "INTT_EventDisplay_new - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    std::cout << "INTT_EventDisplay_new - entry: " << i << std::endl;

    // load entry
    tree->GetEntry(i);
    const auto& clusters = container->clusters();

//     // count clusters
//     const auto tpot_clusters = std::count_if( clusters.begin(), clusters.end(), []( const auto& cluster ){ return cluster._layer == 55 || cluster._layer == 56; } );
//     const auto intt_clusters = std::count_if( clusters.begin(), clusters.end(), []( const auto& cluster ){ return cluster._layer >= 3 && cluster._layer < 7; } );
//
//     if( !(tpot_clusters && intt_clusters) ) continue;

    const auto bco = container->events()[0]._gtm_bco;

    // create canvas
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    cv->cd(1);
    std::unique_ptr<TH2> h0( new TH2I( "h0", "", 200, -95, 95, 200, -95, 95 ) );
    // std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -60, 60, 100, -95, 15 ) );
    h0->SetTitle( Form( "runnumber: %i, event: %i, BCO: %lu", runnumber, i, bco ) );
    h0->GetXaxis()->SetTitle( "x (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->Draw();

    std::unique_ptr<TMarker> marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerColor(2);

    draw_intt_xy();
    draw_tpot_xy();

    // loop over clusters
    for( const auto& cluster:clusters )
    {
      if( cluster._layer != 56 )
      { marker->DrawMarker( cluster._x, cluster._y ); }
    }


    cv->cd(2);
    std::unique_ptr<TH2> h1( new TH2I( "h1", "", 100, -110, 110, 200, -95, 95 ) );
    h1->GetXaxis()->SetTitle( "z (cm)" );
    h1->GetYaxis()->SetTitle( "r (cm)" );
    h1->Draw();

    draw_intt_rz();
    draw_tpot_rz();

    for( const auto& cluster:clusters )
    {

      if( swap_tpc_z && cluster._layer >= 7 && cluster._layer < 55 )
      {
        marker->DrawMarker( -cluster._z,  (cluster._y > 0) ? cluster._r:-cluster._r );
      } else if( cluster._layer != 55 ) {
        marker->DrawMarker( cluster._z,  (cluster._y > 0) ? cluster._r:-cluster._r );
      }
      // { marker->DrawMarker( cluster._z,  cluster._y ); }
    }

    pdfDocument.Add( cv.get() );
  }

}
