#include <RootUtil/PdfDocument.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <optional>
#include <g4eval_hp/MicromegasGeometryContainer.h>

#include <nlohmann/json.hpp>

#include "EvtDisplayUtil.h"

const std::string status = "Internal";


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
  int n_clusters = 0; 
  
  // GlobalX
  std::vector<float>* cluster_x = nullptr;   
  std::vector<float>* cluster_y = nullptr;
  std::vector<float>* cluster_z = nullptr; 
};

template<class T>
  inline constexpr T square( const T& x ) { return x*x; }
    
template<class T>
  inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }


namespace
{
  
  static constexpr double rad_to_deg = 180./M_PI;
  
  template <typename T> 
    inline constexpr int sign(const T& val) 
  { return (T(0) < val) - (val < T(0)); }
  
  // geometry container
  std::unique_ptr<MicromegasGeometryContainer> m_geometry_container;

  // load geometry from file
  void load_geometry( const std::string& filename )
  {
    auto inputfile( TFile::Open( filename.c_str() ) );    
    m_geometry_container.reset( dynamic_cast<MicromegasGeometryContainer*>( inputfile->Get( "MicromegasGeometryContainer" ) ) );
    m_geometry_container->identify();
  }

  // global run number
  int m_runnumber = 0;
  
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
    tree->SetBranchAddress("PhyTrgBCO",&mvtx_event.get()->lvl1_bco);
    tree->SetBranchAddress("Nclusters",&mvtx_event.get()->n_clusters);
    tree->SetBranchAddress("clusterX",&mvtx_event.get()->cluster_x);
    tree->SetBranchAddress("clusterY",&mvtx_event.get()->cluster_y);
    tree->SetBranchAddress("clusterZ",&mvtx_event.get()->cluster_z);
  }
  
  //_______________________________________________
  void setup_mvtx_tree_new( TTree* tree )
  {    
    tree->SetBranchAddress("L1TrigBCO",&mvtx_event.get()->lvl1_bco);
    tree->SetBranchAddress("TotalCleanClus",&mvtx_event.get()->n_clusters);
    tree->SetBranchAddress("CleanClusX",&mvtx_event.get()->cluster_x);
    tree->SetBranchAddress("CleanClusY",&mvtx_event.get()->cluster_y);
    tree->SetBranchAddress("CleanClusZ",&mvtx_event.get()->cluster_z);
  }
  // running tpot event
  FullEventContainer* tpot_event = nullptr;
    
  //_______________________________________________
  void setup_tpot_tree( TTree* tree )
  { tree->SetBranchAddress( "Event", &tpot_event ); }

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
    line->SetLineColorAlpha(kBlack, 0.3);
    const auto INTTstavepointsXY = INTTStavePositionXY();
    for (auto &stave : INTTstavepointsXY)
    { line->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }
  }
  
  //____________________________________________________________________________
  void draw_tpot_xy()
  {
    
    // loop over tiles
    auto line = std::make_unique<TLine>();
    line->SetLineColorAlpha(kBlack, 0.3);
    for( int tile = 0; tile < 8; ++tile )
    {
      const auto first_strip_center = get_strip_center( 55, tile, 0 );
      const auto last_strip_center = get_strip_center( 55, tile, 255 );
      line->DrawLine( first_strip_center.x(), first_strip_center.y(), last_strip_center.x(), last_strip_center.y() );
    }
  }
  
  //____________________________________________________________________________
  void draw_mvtx_xy()
  {  
    auto line = std::make_unique<TLine>();
    line->SetLineColorAlpha(kBlack, 0.3);
    const auto MVTXstavepointsXY = MVTXStavePositionXY();
    for (auto &stave : MVTXstavepointsXY)
    { line->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }
  }
  
  //____________________________________________________________________________
  void draw_intt_rz()
  {
    const auto INTTstavepoints_ZRho_min = INTTStaveMinPositionZRho();
    const auto INTTstavepoints_ZRho_max = INTTStaveMaxPositionZRho();
    auto line_min = std::make_unique<TLine>();
    line_min->SetLineColorAlpha(kBlack, 0.3);
    
    auto line_max = std::make_unique<TLine>();
    line_max->SetLineColorAlpha(kBlack, 0.3);
    line_max->SetLineStyle(3);
    for(const auto &stave : INTTstavepoints_ZRho_min)
    { line_min->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }

    for(const auto &stave : INTTstavepoints_ZRho_max)
    { line_max->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }
  }
  
  //____________________________________________________________________________
  void draw_tpot_rz()
  {
    
    // loop over tiles
    auto line = std::make_unique<TLine>();
    line->SetLineColorAlpha(kBlack, 0.3);
    for( int tile = 0; tile < 8; ++tile )
    {
      const auto first_strip_center = get_strip_center( 56, tile, 0 );
      const auto last_strip_center = get_strip_center( 56, tile, 255 );
      line->DrawLine( 
        first_strip_center.z(), -get_r(first_strip_center.x(), first_strip_center.y()),
        last_strip_center.z(), -get_r(last_strip_center.x(), last_strip_center.y()) );
    }
  }
  
  //____________________________________________________________________________
  void draw_mvtx_rz()
  {
    const auto stavepoints_ZRho_min = MVTXStaveMinPositionZRho();
    const auto stavepoints_ZRho_max = MVTXStaveMaxPositionZRho();
    auto line_min = std::make_unique<TLine>();
    line_min->SetLineColorAlpha(kBlack, 0.3);
    
    auto line_max = std::make_unique<TLine>();
    line_max->SetLineColorAlpha(kBlack, 0.3);
    line_max->SetLineStyle(3);
    
    for( const auto &stave : stavepoints_ZRho_min)
    { line_min->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }
    
    for( const auto &stave : stavepoints_ZRho_max)
    { line_max->DrawLine(get<0>(stave), get<1>(stave), get<2>(stave), get<3>(stave)); }
  }
  
  // track parameter
  class track_parameters_t
  {
    public:
    
    double xv = 0;
    double yv = 0;
    double zv = 0;
    
    double px = 0;
    double py = 0;
    double pz = 0;
    
    double length = 0;
    int q = 1;
    
    using list = std::vector<track_parameters_t>;
  
  };
  
  track_parameters_t::list tracks;
  
  //____________________________________________________________________________
  void write_to_json( const TString& filename )
  {
    
    std::cout << "write_to_json - filename: " << filename << std::endl;

    nlohmann::json runstats;
    runstats.push_back( "sPHENIX Experiment at RHIC" );
    runstats.push_back( "Data recorded: 2023-08-11" );
    runstats.push_back( Form( "Run / Event: %i/%i", m_runnumber,  intt_event->event_id ) );
    runstats.push_back( Form( "BCO: %lld",  intt_event->lvl1_bco ) );
    runstats.push_back( "Cosmics" );
    nlohmann::json event;
    event["runstats"] = runstats;
    event["B"] = 1.4;
    
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
    
    // add mvtx clusters
    for( int iclus = 0; iclus < mvtx_event->n_clusters; ++iclus )
    { 
      nlohmann::json cluster = {
        {"x", (*mvtx_event->cluster_x)[iclus]},
        {"y", (*mvtx_event->cluster_y)[iclus]},
        {"z", (*mvtx_event->cluster_z)[iclus]},
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
   
    // add tracks
    nlohmann::json jstracks;    
    jstracks["B"]=0.00001;
    for( const auto& track:tracks )
    {
      nlohmann::json jstrack =  
      {
        { "color", 16776960 },
        { "xyz", { track.xv, track.yv, track.zv } },
        { "pxyz", { track.px, track.py, track.pz } },
        { "l", track.length },
        { "nh", 50 },
        {"q", track.q }
      };
    
      jstracks["TRACKHITS"].emplace_back( jstrack );
    }
      
    nlohmann::json main;
//     main["EVENT"]=event;
//     main["META"]=meta;
//     main["HITS"]=hits;            
    
    if( !tracks.empty() ) { main["TRACKS"]=jstracks; }
    
    // write event to file
    std::ofstream out( filename.Data() );
    out << std::setw(2) << main << std::endl;
    out.close();  
        
  }

  
  //_________________________________________________________________________________
  using position_t = std::array<double,2>;
  using position_vector_t = std::vector<position_t>;
  
  using line_fit_output_t = std::array<double,2>;
  line_fit_output_t line_fit(const position_vector_t& positions)
  {
    double xsum=0;
    double x2sum=0;
    double ysum=0;
    double xysum=0;
    for( const auto& [r,z]:positions )
    {
      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+square(r);              //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
    
    const auto npts = positions.size();
    const double denominator = (x2sum*npts-square(xsum));
    const double a= (xysum*npts-xsum*ysum)/denominator;            //calculate slope
    const double b= (x2sum*ysum-xsum*xysum)/denominator;           //calculate intercept
    return { a, b };
  }

  using mask_list_t = std::set<size_t>;
  using mask_map_t = std::map<Long64_t, mask_list_t>;
  
  mask_map_t hit_mask;
  mask_map_t fit_hit_mask = 
  {
    { 126530705406, {1} },
    { 128292557799, {2} },
    { 128434167340, {6, 7} }, 
    { 128755275235, {0, 1}} 
  };
  
  //______________________________________________________________________
  position_vector_t apply_mask(Long64_t lvl1_bco, const position_vector_t& source, const mask_map_t& mask_map  )
  {
    const auto iter = mask_map.find( lvl1_bco );
    if( iter == mask_map.end() ) return source;
    const auto& mask = iter->second;
    
    position_vector_t out;
    
    for( size_t i = 0; i < source.size(); ++i )
    {  
      if( mask.find( i ) == mask.end() ) out.push_back( source[i] ); 
    
    }

    return out;
    
  }
    
  //______________________________________________________________________
  position_vector_t get_positions_xy(Long64_t lvl1_bco)
  {   
  
    position_vector_t positions;

    // INTT clusters
    for( int iclus = 0; iclus < intt_event->n_clusters; ++iclus )
    { positions.push_back( { (*intt_event->cluster_x)[iclus]/10, (*intt_event->cluster_y)[iclus]/10} ); }
    
    // MVTX hits
    for( size_t iclus = 0; iclus < mvtx_event->n_clusters; ++iclus )
    { 
      // mask noisy mvtx cluster manually
      if( std::abs((*mvtx_event->cluster_x)[iclus]+4.05874)<0.01 && std::abs((*mvtx_event->cluster_y)[iclus]-0.300536)<0.01 )
      { continue; }
      
      positions.push_back( { (*mvtx_event->cluster_x)[iclus], (*mvtx_event->cluster_y)[iclus] } ); 
    }
    
    // TPOT clusters
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 55 )
      { positions.push_back( { cluster.center.x(), cluster.center.y() } ); }
    } 
    
    std::cout << "get_positions_xy - bco: " << lvl1_bco << std::endl;
    int index = 0;
    for( const auto& [x,y]:positions )
    { std::cout << " index: " << index++ << " (" << x << ", " << y << ")" << std::endl; } 
    
    return positions;    
    
  }
  
  //______________________________________________________________________
  position_vector_t get_positions_rz(Long64_t lvl1_bco)
  {
   
    position_vector_t positions;

    // INTT clusters
    for( int iclus = 0; iclus < intt_event->n_clusters; ++iclus )
    { 
      positions.push_back( 
      {
        (*intt_event->cluster_r)[iclus]/10, 
        (*intt_event->cluster_z)[iclus]/10
      } ); 
    }
    
    // MVTX hits
    for( size_t iclus = 0; iclus < mvtx_event->n_clusters; ++iclus )
    { 
      // mask noisy mvtx cluster manually
      if( std::abs((*mvtx_event->cluster_x)[iclus]+4.05874)<0.01 && std::abs((*mvtx_event->cluster_y)[iclus]-0.300536)<0.01 )
      { continue; }

      positions.push_back(
      {
        sign( (*mvtx_event->cluster_y)[iclus] )*get_r( (*mvtx_event->cluster_x)[iclus], (*mvtx_event->cluster_y)[iclus] ),
        (*mvtx_event->cluster_z)[iclus]
      } );
    }
    
    // TPOT clusters
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 56 )
      { 
        positions.push_back( 
        {
          -get_r(cluster.center.x(), cluster.center.y()),
          cluster.center.z() 
        } ); 
      }
    }
    
    
    std::cout << "get_positions_rz - bco: " << lvl1_bco << std::endl;
    int index = 0;
    for( const auto& [r,z]:positions )
    { std::cout << " index: " << index++ << " (" << r << ", " << z << ")" << std::endl; } 

    return positions;    
  }
  
  template< class T> 
    void print_list( const std::string& name, const std::string& type, const std::vector<T>& list )
  {
    
    std::cout << "const std::vector<" <<  type << "> " << name << " = {" << std::endl;
    bool first = true;
    int count = 0;
    for( const auto& value:list )
    {
      if( !first ) std::cout << ", ";
      first = false;
      if( count == 10 ) 
      {
        count = 0;
        std::cout << std::endl;
      }
      std::cout << " " << value;
      ++count;
    }
    std::cout << std::endl << "};" << std::endl;
  }
  
}

//____________________________________________________________________________
TString MVTX_EventDisplay_nofield( int runnumber = 25926 )
{
  
  // store run number globally. Needed for JSon
  m_runnumber = runnumber;
  
  // load geometry
  load_geometry( "micromegas_geometry.root" );

  // pdf output
  const TString pdffilename = Form( "Figures/MVTX_event_display-%08i-0000.pdf", runnumber );
  PdfDocument pdfDocument( pdffilename );
  
  // intt input file
  const TString intt_inputfilename = Form( "/sphenix/user/ChengWei/INTT/INTT_commissioning/cosmic/%i/INTT_eventdisplay_cluster_%i.root", runnumber, runnumber );
  auto intt_tfile = TFile::Open( intt_inputfilename, "READ" );
  auto intt_tree = static_cast<TTree*>( intt_tfile->Get("tree_clu") );
  setup_intt_tree( intt_tree );
  
  // mvtx input file
  // const TString mvtx_inputfilename = "/sphenix/user/zshi/NewCameron/MVTXClus/MVTXClusFile.root";
  const TString mvtx_inputfilename = "/sphenix/user/zshi/NewCameron/MVTXClus/Run25926.root";
  auto mvtx_tfile = TFile::Open( mvtx_inputfilename, "READ" );
  auto mvtx_tree = static_cast<TTree*>( mvtx_tfile->Get("MVTXClusTree") );
  setup_mvtx_tree_new( mvtx_tree );
    
  // tpot input file
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-full_filtered.root", runnumber );
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  setup_tpot_tree( tpot_tree );

  const auto intt_entries = intt_tree->GetEntries();
  // const auto intt_entries = 1000;
  const auto mvtx_entries = mvtx_tree->GetEntries();
  const auto tpot_entries = tpot_tree->GetEntries();
  
  std::cout << "MVTX_EventDisplay - intt_inputfilename: " << intt_inputfilename << " entries: " << intt_entries << std::endl;
  std::cout << "MVTX_EventDisplay - mvtx_inputfilename: " << mvtx_inputfilename << " entries: " << mvtx_entries<< std::endl;
  std::cout << "MVTX_EventDisplay - tpot_inputfilename: " << tpot_inputfilename << " entries: " << tpot_entries << std::endl;
  std::cout << "MVTX_EventDisplay - pdffilename: " << pdffilename << std::endl;
  
  int intt_entry = 0;
  int mvtx_entry = 0;
  int tpot_entry = 0;
  
  std::vector<Long64_t> accepted_bco;

  std::set<Long64_t> lvl1_mask;
  // std::set<Long64_t> lvl1_mask = { 881044278792 };
  // std::set<Long64_t> lvl1_mask = { 881044278792,  881792317048, 883479809083, 883639102957,  885501441964,  885740684584,  887523516007, 888261908398 };  
  for(; intt_entry < intt_entries && mvtx_entry < mvtx_entries && tpot_entry < tpot_entries; ++intt_entry )
  {
    // load intt_entry
    intt_tree->GetEntry(intt_entry);
    
    // skip if no cluster
    if( intt_event->n_clusters <= 0 ) continue;
    
    // get bco
    const auto& intt_bco = intt_event->lvl1_bco;
    
    if( !lvl1_mask.empty() && lvl1_mask.find( intt_bco ) == lvl1_mask.end() ) continue;
    
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

    auto positions_xy = get_positions_xy(intt_bco);
    auto masked_positions_xy = apply_mask( intt_bco, positions_xy, hit_mask );
    auto fit_positions_xy = apply_mask( intt_bco, positions_xy, fit_hit_mask );
    std::sort( positions_xy.begin(), positions_xy.end(), []( const position_t& first, const position_t& second ) { return first[1] < second[1]; } );
    std::sort( fit_positions_xy.begin(), fit_positions_xy.end(), []( const position_t& first, const position_t& second ) { return first[1] < second[1]; } );
    
    auto positions_rz = get_positions_rz(intt_bco);
    auto masked_positions_rz = apply_mask( intt_bco, positions_rz, hit_mask );
    auto fit_positions_rz = apply_mask( intt_bco, positions_rz, fit_hit_mask );
    std::sort( positions_rz.begin(), positions_rz.end(), []( const position_t& first, const position_t& second ) { return first[0] < second[0]; } );
    std::sort( fit_positions_rz.begin(), fit_positions_rz.end(), []( const position_t& first, const position_t& second ) { return first[0] < second[0]; } );
       
    // create canvas
    // std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 1960, 900 ) );
    auto cv( new TCanvas( "cv", "cv", 1960, 900 ) );
    cv->Divide( 2, 1 );

    cv->cd(1);
    gPad->SetTopMargin( 0.07 );

    // std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -15, 15, 100, -15, 15 ) );
    auto h0( new TH2I( "h0", "", 100, -60, 60, 100, -95, 15 ) );
    // h0->SetTitle( Form( "runnumber: %i, event: %i, BCO: %lld", runnumber, intt_event->event_id, intt_bco ) );
    h0->SetTitle( "" );
    h0->GetXaxis()->SetTitle( "x [cm]" );
    h0->GetYaxis()->SetTitle( "y [cm]" );
    h0->Draw();
    
    draw_intt_xy();
    draw_tpot_xy();
    draw_mvtx_xy();

    {
      auto marker( new TMarker() );
      marker->SetMarkerStyle(20);
      marker->SetMarkerColor(1);    
      for( const auto& [x,y]:masked_positions_xy )
      { marker->DrawMarker( x, y ); }
    }
    
    {
      auto marker( new TMarker() );
      marker->SetMarkerStyle(20);
      marker->SetMarkerColor(2);    
      for( const auto& [x,y]:fit_positions_xy )
      { marker->DrawMarker( x, y ); }
    }
    
    // perform linear fit
    const auto [slope_xy,y0] = line_fit( fit_positions_xy );
    {
      const auto [xf,yf] = fit_positions_xy.front();
      const auto xbegin = (yf-y0)/slope_xy;
      // const auto ybegin = y0 + slope_xy*xf;
      
      const auto [xb,yb] = fit_positions_xy.back();
      const auto xend = (yb-y0)/slope_xy;
      // const auto yend = y0 + slope_xy*xb;
      
      auto line = new TLine( xbegin, yf, xend, yb );
      // auto line = new TLine( xf, ybegin, xb, yend );
      line->SetLineColor(4);
      line->Draw();
    }
    
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      // text->DrawLatex( 0.77, 0.95, "#it{08/11/2023}" );
      text->DrawLatex( 0.12, 0.95, Form( "#it{08/23/2023} run: %i, BCO: %lld", runnumber, intt_bco ) );
      // text->DrawLatex( 0.1, 0.95, Form( "run: %i, event: %i, BCO: %lld", runnumber, intt_event->event_id, intt_bco ) );
    }
    
    {
      auto text = new TPaveText(0.66,0.74,0.94,0.89, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      text->AddText("Cosmic data");
      text->AddText(Form( "Run %i", runnumber ) );
            
      text->Draw();
    }


    cv->cd(2);
    gPad->SetTopMargin( 0.07 );
    auto h1( new TH2I( "h1", "", 100, -110, 110, 100, -95, 15 ) );
    h1->GetXaxis()->SetTitle( "z [cm]" );
    h1->GetYaxis()->SetTitle( "r [cm]" );
    h1->Draw();
    
    draw_intt_rz();
    draw_tpot_rz();
    draw_mvtx_rz();

    {
      auto marker( new TMarker() );
      marker->SetMarkerStyle(20);
      marker->SetMarkerColor(1);    
      for( const auto& [r,z]:masked_positions_rz )
      { marker->DrawMarker(z,r); }
    }
    
    {
      auto marker( new TMarker() );
      marker->SetMarkerStyle(20);
      marker->SetMarkerColor(2);    
      for( const auto& [r,z]:fit_positions_rz )
      { marker->DrawMarker(z,r); }
    }
    
    const auto [slope,z0] = line_fit( fit_positions_rz );
    {
      const auto [rf,zf] = fit_positions_rz.front();
      const auto zbegin = z0 + slope*rf;
      
      const auto [rb,zb] = fit_positions_rz.back();
      const auto zend = z0 + slope*rb;
      
      auto line = new TLine( zbegin, rf, zend, rb );
      line->SetLineColor(4);
      line->Draw();
    }
    
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      // text->DrawLatex( 0.77, 0.95, "#it{08/11/2023}" );
      text->DrawLatex( 0.12, 0.95, Form( "#it{08/23/2023} run: %i, BCO: %lld", runnumber, intt_bco ) );
      // text->DrawLatex( 0.1, 0.95, Form( "run: %i, event: %i, BCO: %lld", runnumber, intt_event->event_id, intt_bco ) );
    }
    
    {
      auto text = new TPaveText(0.66,0.74,0.94,0.89, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      text->AddText("Cosmic data");
      text->AddText(Form( "Run %i", runnumber ) );
      text->Draw();
    }

    pdfDocument.Add( cv );
    const TString png_filename = Form( "Figures/MVTX_event_display-%08i-%lld.png", runnumber, intt_bco );
    const TString pdf_filename = Form( "Figures/MVTX_event_display-%08i-%lld.pdf", runnumber, intt_bco );
    cv->SaveAs( png_filename );
    cv->SaveAs( pdf_filename );
   
    // clear tracks
    tracks.clear();
    
    if( true )
    {
      const auto [xb,yb] = fit_positions_xy.back();
      const double rb = get_r(xb, yb);
      const double zb = z0+slope*rb;
      
      const auto [xf,yf] = fit_positions_xy.front();
      const double rf = get_r(xf, yf);
      const double zf = z0 - slope*rf;      
      
      const auto phi_xy = std::atan( slope_xy );
      const double pt = 50;

      track_parameters_t track;
            
      // get vertex position
      track.xv = xf;
      track.yv = yf;
      track.zv = zf;

      // track momentum
      track.px = pt*std::cos( phi_xy );
      track.py = pt*std::sin( phi_xy );
      track.pz = pt*slope;
      
      // make sure that py is negative
      if( track.py < 0 ) 
      {
        track.px *= -1;
        track.py *= -1;
      }

      // track.length = 2.*std::sqrt( square( xf-xb ) + square( yf-yb ) + square( zf-zb ) );
      track.length = 170.;
      tracks.push_back( track );
    }
    
    // save
    accepted_bco.push_back( intt_bco );
    
    // generate json file
    const auto json_filename = Form( "event_display/mvtx_event_display_%08i_%08i_%lld.json", runnumber, intt_event->event_id, intt_bco );
    write_to_json( json_filename);
  }

  print_list( "accepted_bco", "uint64_t", accepted_bco );
  
  return pdffilename;
  
}
