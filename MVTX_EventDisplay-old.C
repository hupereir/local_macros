#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <optional>
#include <g4eval_hp/MicromegasGeometryContainer.h>

#include <nlohmann/json.hpp>

// TPOT cluster
class Cluster
{
  public:

  unsigned short layer = 0;
  
  unsigned short tile = 0;

  unsigned short size = 0;
  
  double charge = 0;

  int strip = 0;

  int region = 0;
  
  TVector3 begin;
  TVector3 end;  
  TVector3 center;  

  bool is_signal = true;
  
  using List = std::vector<Cluster>;
};

// TPOT event
class FullEventContainer: public TObject
{
  
  public:
  
  /// ll1 bco
  uint64_t lvl1_bco = 0;

  /// ll1 counter
  uint32_t lvl1_counter = 0;
  
  // number of waveforms
  unsigned short n_waveforms_all = 0;
  
  // number of signal waveforms
  unsigned short n_waveforms_signal = 0;
  
  // number of clusters
  unsigned short n_clusters = 0;
  
  // clusters
  Cluster::List clusters;
  
  // number of clusters per detector
  std::vector<unsigned short> n_detector_clusters;
 
  // number of clusters per region
  std::vector<unsigned short> n_region_clusters;
 
  // minimum cluster charge per detector
  std::vector<unsigned short> min_cluster_size;
 
  // minimum cluster charge per detector
  std::vector<double> min_cluster_charge;
  
  // number of clusters
  unsigned short n_phi_clusters = 0;
  
  // number of clusters
  unsigned short n_z_clusters = 0;
    
  ClassDef(FullEventContainer,1)
  
};

//! lines
using line_t = std::array<double,4>;

class InttEvent
{
  public:
  
  int event_id = 0;
  Long64_t lvl1_bco = 0;
  int n_clusters = 0;
  std::vector<double>* cluster_x = nullptr;
  std::vector<double>* cluster_y = nullptr;
  std::vector<double>* cluster_z = nullptr;
  std::vector<double>* cluster_r = nullptr;
};


class MvtxEvent
{
  public: 
  
  // L1TrigBCOSave
  int64_t lvl1_bco = 0;
  
  // NumberHits
  int n_hits = 0; 
  
  // GlobalX
  std::vector<float>* cluster_x = nullptr;   
  std::vector<float>* cluster_y = nullptr;
  std::vector<float>* cluster_z = nullptr; 
};

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

  // runnint intt event
  auto intt_event = std::make_unique<InttEvent>();
  
  //_______________________________________________
  void setup_intt_tree( TTree* tree )
  {    
    tree->SetBranchAddress("eID",&intt_event.get()->event_id);
    tree->SetBranchAddress("bco_full",&intt_event.get()->lvl1_bco);
    tree->SetBranchAddress("nclu",&intt_event.get()->n_clusters);
    tree->SetBranchAddress("clu_x",&intt_event.get()->cluster_x);
    tree->SetBranchAddress("clu_y",&intt_event.get()->cluster_y);
    tree->SetBranchAddress("clu_z",&intt_event.get()->cluster_z);
    tree->SetBranchAddress("clu_r",&intt_event.get()->cluster_r);
  }
  
  // running mvtx event
  auto mvtx_event = std::make_unique<MvtxEvent>();
  
  //_______________________________________________
  void setup_mvtx_tree( TTree* tree )
  {    
    tree->SetBranchAddress("L1TrigBCOSave",&mvtx_event.get()->lvl1_bco);
    tree->SetBranchAddress("NumberHits",&mvtx_event.get()->n_hits);
    tree->SetBranchAddress("GlobalX",&mvtx_event.get()->cluster_x);
    tree->SetBranchAddress("GlobalY",&mvtx_event.get()->cluster_y);
    tree->SetBranchAddress("GlobalZ",&mvtx_event.get()->cluster_z);
  }
  
  // running tpot event
  FullEventContainer* tpot_event = nullptr;
    
  //_______________________________________________
  void setup_tpot_tree( TTree* tree )
  { tree->SetBranchAddress( "Event", &tpot_event ); }

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
  
  
  //_______________________________________________
  TVector3 get_strip_center( int layer, int tile, int strip )
  {
    // update cluster position
    const auto strip_begin = m_geometry_container->get_strip_begin( layer, tile, strip );
    const auto strip_end = m_geometry_container->get_strip_end( layer, tile, strip );
    return 0.5*(strip_begin + strip_end);
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
  void draw_intt_rz()
  {
    for( const double& r:{-7.7, -10.3, 7.7, 10.3} )
    {
      auto line = new TLine( -25, r, 25, r );
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

  //____________________________________________________________________________
  void write_to_json( const TString& filename )
  {

    // prepare json header    
    // event
    nlohmann::json event = 
    {
      { "runid", 1 },
      { "evtid", intt_event->event_id },
      { "time", 0 },
      { "type", "Cosmics" },
      { "s_nn", 0 },
      { "B", 1.4 },
      { "pv", {0,0,0}},
    };
    
    // meta data
    nlohmann::json meta = 
    { { "HITS",{
        { "INNERTRACKER", 
          {
            { "type", "3D" },
            { "options", { {"size", 10}, { "color", 16777215 } } }
          }
        },
        { "TRACKHITS",
          {
            { "type", "3D" },
            { "options", { {"size", 2 }, { "transparent", 0.8 }, { "color", 16777215 } } }
          }
        }
      } }
    };
        
    nlohmann::json hits;    
    
    // add intt clusters
    for( int iclus = 0; iclus < intt_event->n_clusters; ++iclus )
    { 
      nlohmann::json cluster = {
        {"x", (*intt_event->cluster_x)[iclus]/10},
        {"y", (*intt_event->cluster_y)[iclus]/10},
        {"z", (*intt_event->cluster_z)[iclus]/10},
        {"e", 1 }
      };
      hits["INNERTRACKER"].emplace_back( cluster );
    }
    
    // add TPOT clusters
    static constexpr int npoints = 100;
    for( const auto& cluster:tpot_event->clusters )
    {
      for( int ipoint = 0; ipoint <= npoints; ++ipoint )
      {
        double fraction = double(ipoint)/npoints;
        const auto strip_center = fraction*cluster.begin + (1.0-fraction)*cluster.end;
     
        nlohmann::json strip_info = {
          {"x", strip_center.x()},
          {"y", strip_center.y()},
          {"z", strip_center.z()},
          {"e", 1 }
        };
        hits["TRACKHITS"].emplace_back( strip_info );
      }
      
    }
   
    nlohmann::json main;
    main["EVENT"]=event;
    main["META"]=meta;
    main["HITS"]=hits;            
  
    // write event to file
    std::ofstream out( filename.Data() );
    out << std::setw(2) << main << std::endl;
    out.close();  
        
  }
}

//____________________________________________________________________________
TString MVTX_EventDisplay( int runnumber = 25475 )
{
  
  // load intt geometry
  load_intt_geometry( "INTT_TLine_geo.root" );
  
  // load geometry
  load_geometry( "micromegas_geometry.root" );

  // pdf output
  const TString pdffilename = Form( "Figures/MVTX_event_display-%08i-0000.pdf", runnumber );
  PdfDocument pdfDocument( pdffilename );
  
  // intt input file
  // const TString intt_inputfilename = Form( "/sphenix/user/ChengWei/INTT/INTT_commissioning/cosmic/%i/INTT_eventdisplay_cluster.root", runnumber );
  const TString intt_inputfilename = Form( "/sphenix/user/ChengWei/INTT/INTT_commissioning/cosmic/%i/new_survey_geo/INTT_eventdisplay_cluster.root", runnumber );
  auto intt_tfile = TFile::Open( intt_inputfilename, "READ" );
  auto intt_tree = static_cast<TTree*>( intt_tfile->Get("tree_clu") );
  setup_intt_tree( intt_tree );
  
  // mvtx input file
  // const TString mvtx_inputfilename = Form( "/sphenix/user/zshi/TestBuild/FastMLUConn/MVTXINTTEvtDisplay/macros/detectors/sPHENIX/OutFile/Run%i/3DHitTree_0.root", runnumber );
  const TString mvtx_inputfilename = Form( "/sphenix/user/hjheng/macros/detectors/sPHENIX/OutFile/ClusterTree/ClusterTree_Run%i_0.root", runnumber );

  auto mvtx_tfile = TFile::Open( mvtx_inputfilename, "READ" );
  auto mvtx_tree = static_cast<TTree*>( mvtx_tfile->Get("MVTXPixelTree") );
  setup_mvtx_tree( mvtx_tree );
  
  // tpot input file
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000.root", runnumber );
  // const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  setup_tpot_tree( tpot_tree );
  
  std::cout << "MVTX_EventDisplay - intt_inputfilename: " << intt_inputfilename << std::endl;
  std::cout << "MVTX_EventDisplay - mvtx_inputfilename: " << mvtx_inputfilename << std::endl;
  std::cout << "MVTX_EventDisplay - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  std::cout << "MVTX_EventDisplay - pdffilename: " << pdffilename << std::endl;

  const auto intt_entries = intt_tree->GetEntries();
  // const auto intt_entries = 1000;
  const auto mvtx_entries = mvtx_tree->GetEntries();
  const auto tpot_entries = tpot_tree->GetEntries();
  
  int intt_entry = 0;
  int mvtx_entry = 0;
  int tpot_entry = 0;
  for(; intt_entry < intt_entries && mvtx_entry < mvtx_entries && tpot_entry < tpot_entries; ++intt_entry )
  {
    // load intt_entry
    intt_tree->GetEntry(intt_entry);
    
    // skip if no cluster
    if( intt_event->n_clusters <= 0 ) continue;
    
    // get bco
    const auto& intt_bco = intt_event->lvl1_bco;
    
    // find matching tpot entry
    uint64_t tpot_bco = 0;
    while( tpot_entry < tpot_entries )
    {
      tpot_tree->GetEntry( tpot_entry );
      tpot_bco = tpot_event->lvl1_bco;
     
      if( tpot_bco >= intt_bco ) break;
      
      ++tpot_entry;
    }
        
    // check that both BCO match within 2 units
    if( std::abs<int64_t>( tpot_bco - intt_bco ) > 1 )
    {
      std::cout << "MVTX_EventDisplay -"
        << " intt_entry: " << intt_entry 
        << " intt_bco: " << intt_bco
        << " not found (TPOT)"
        << std::endl;
      continue;
    }
    
    // find matching mvtx entry
    int64_t mvtx_bco = 0;
    while( mvtx_entry < mvtx_entries )
    {
      mvtx_tree->GetEntry( mvtx_entry );
      mvtx_bco = mvtx_event->lvl1_bco;
     
      if( mvtx_bco >= intt_bco ) break;
      
      ++mvtx_entry;
    }
        
    // check that both BCO match within 2 units
    if( std::abs<int64_t>( mvtx_bco - intt_bco ) > 1 )
    {
      std::cout << "MVTX_EventDisplay -"
        << " intt_entry: " << intt_entry 
        << " intt_bco: " << intt_bco
        << " not found (MVTX)"
        << std::endl;
      continue;
    }
    
    
    std::cout << "MVTX_EventDisplay -"
      << " intt_entry: " << intt_entry 
      << " intt_bco: " << intt_bco
      << " tpot_bco: " << tpot_bco
      << " mvtx_bco: " << mvtx_bco
      << std::endl;

    // create canvas
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    cv->cd(1);

    // std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -15, 15, 100, -15, 15 ) );
    std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -60, 60, 100, -95, 15 ) );
    h0->SetTitle( Form( "runnumber: %i, event: %i, BCO: %lld", runnumber, intt_event->event_id, intt_bco ) );
    h0->GetXaxis()->SetTitle( "x (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->Draw();
    
    std::unique_ptr<TMarker> marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerColor(2);
    
    // draw INTT clusters
    draw_intt_xy();
    for( int iclus = 0; iclus < intt_event->n_clusters; ++iclus )
    { marker->DrawMarker( (*intt_event->cluster_x)[iclus]/10, (*intt_event->cluster_y)[iclus]/10 ); }
    
    // draw MVTX hits
    const size_t n_hits = (*mvtx_event->cluster_x).size();
    std:: cout << "MVTX_EventDisplay n_hits: " << n_hits << std::endl;
    for( size_t ihit = 0; ihit < n_hits; ++ihit )
    { 
      std:: cout << "MVTX_EventDisplay cluster_x: " 
        << (*mvtx_event->cluster_x)[ihit] << ", "
        << (*mvtx_event->cluster_y)[ihit] 
        << std::endl;    
      marker->DrawMarker( (*mvtx_event->cluster_x)[ihit], (*mvtx_event->cluster_y)[ihit] ); 
    }
    
    // draw TPOT clusters
    draw_tpot_xy();
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 55 )
      { marker->DrawMarker( cluster.center.x(), cluster.center.y() ); }
    }
        
    cv->cd(2);
    // std::unique_ptr<TH2> h1( new TH2I( "h1", "", 100, -50, 50, 100, -15, 15 ) );
    std::unique_ptr<TH2> h1( new TH2I( "h1", "", 100, -110, 110, 100, -95, 15 ) );
    h1->GetXaxis()->SetTitle( "z (cm)" );
    h1->GetYaxis()->SetTitle( "r (cm)" );
    h1->Draw();
    
    draw_intt_rz();
    for( int iclus = 0; iclus < intt_event->n_clusters; ++iclus )
    { marker->DrawMarker( (*intt_event->cluster_z)[iclus]/10, (*intt_event->cluster_r)[iclus]/10 ); }
        
    // draw TPOT clusters  
    draw_tpot_rz();
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 56 )
      { marker->DrawMarker( cluster.center.z(), -get_r(cluster.center.x(), cluster.center.y()) ); }
    }

    pdfDocument.Add( cv.get() );
    
    // generate json file
    const auto json_filename = Form( "event_display/tpot_intt_event_display_%08i_%08i.json", runnumber, intt_event->event_id );
    write_to_json( json_filename);
  }
  
  return pdffilename;
  
}
