#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>

#include <memory>
#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasRawDataEvaluation.h>
#include <g4eval_hp/MicromegasGeometryContainer.h>


// true if one wants to filter events with at least one cluster in both views of the same module
/* usefull to select cosmic tracks */
bool m_apply_event_filter = false;  

//____________________________________________________________________________
bool operator < (const MicromegasRawDataEvaluation::Waveform& lhs, const MicromegasRawDataEvaluation::Waveform& rhs )
{ return lhs.strip < rhs.strip; }

//____________________________________________________________________________
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

//____________________________________________________________________________
class FullEventContainer: public TObject
{
  
  public:
  
  void Reset()
  {
    lvl1_bco = 0;
    lvl1_counter = 0;
    n_waveforms_all = 0;
    n_waveforms_signal = 0;
    n_clusters = 0;
    n_phi_clusters = 0;
    n_z_clusters = 0;
    clusters.clear();
    n_detector_clusters.clear();
    n_region_clusters.clear();
    min_cluster_size.clear();
    min_cluster_charge.clear();
    first_cluster_strip.clear();
  }
  
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
  
  // strip of the first cludster, per detector
  std::vector<int> first_cluster_strip;
  
  // number of clusters
  unsigned short n_phi_clusters = 0;
  
  // number of clusters
  unsigned short n_z_clusters = 0;
    
  ClassDef(FullEventContainer,1)
  
};

namespace
{

  // counters
  int m_bco_count = 0;
  int m_bco_count_filtered = 0;
  
  // filtered bco
  std::vector<uint64_t> m_filtered_lvl1_bco;
    
  using waveform_set_t = std::set<MicromegasRawDataEvaluation::Waveform>;
  using waveform_map_t = std::map<int,waveform_set_t>;
  
  // number of sigma above pedestal to remove noise
  static constexpr double n_sigma = 5.;
  // static constexpr double n_sigma = 4.;
  // static constexpr double n_sigma = 4.5;
   
  // geometry container
  std::unique_ptr<MicromegasGeometryContainer> m_geometry_container;

  // load geometry from file
  void load_geometry( const std::string& filename )
  {
    auto inputfile( TFile::Open( filename.c_str() ) );    
    m_geometry_container.reset( dynamic_cast<MicromegasGeometryContainer*>( inputfile->Get( "MicromegasGeometryContainer" ) ) );
    m_geometry_container->identify();
  }
  
  // output tfile
  std::unique_ptr<TFile> tfile_out;
 
  // detector names
  // running container
  FullEventContainer* m_fullEventContainer = nullptr;

  // output tree
  TTree* tree_out = nullptr;
  
  // calibration data
  MicromegasCalibrationData calibration_data;
   
  // create output tree
  void create_tree( const TString& rootfilename )
  {
    // create output TFile
    tfile_out.reset(TFile::Open(rootfilename, "RECREATE"));
    tree_out = new TTree( "T", "T" );
    m_fullEventContainer = new FullEventContainer;
    tree_out->Branch( "Event", &m_fullEventContainer );
  }
  
  // save output tree
  void save_tree()
  {
    tfile_out->cd();
    tree_out->Write();
    tfile_out->Close();
  }
  
  // sample windows to determine signal and background clusters
  using sample_window_t = std::pair<unsigned short, unsigned short>;
  sample_window_t sample_window_signal = {20, 45};
  // sample_window_t sample_window_signal = {15, 40};
  sample_window_t sample_window_background = {5, 15 };
      
  //_________________________________________________________________
  Cluster::List get_clusters( const waveform_set_t& waveform_set )
  {
    
    using cluster_t = std::vector<MicromegasRawDataEvaluation::Waveform>;
    using cluster_list_t = std::vector<cluster_t>;
   
    // rolling cluster
    cluster_list_t clusters;
    cluster_t cluster;
    int previous_strip = -1;
    
    for( const auto& waveform:waveform_set )
    {
      if( previous_strip >= 0 && waveform.strip > previous_strip+1 )
      { 
        if( !cluster.empty() ) clusters.push_back(cluster);
        cluster.clear();
      }
      
      previous_strip = waveform.strip;
      cluster.push_back( waveform );
    }
    
    // store last cluster
    if( !cluster.empty() ) clusters.push_back( cluster );
    
    Cluster::List out;
    
    // assign strip and region to cluster
    for( auto&& wf_list:clusters )
    {
      
      assert( !sw_list.empty() );
      
      Cluster cluster;
      cluster.layer = wf_list.front().layer;
      cluster.tile = wf_list.front().tile;
      cluster.size = wf_list.size();
      
      int adc_max = -1;
      for( const auto& waveform:wf_list) 
      {
        if( waveform.adc_max > adc_max ) 
        { 
          cluster.strip = waveform.strip; 
          adc_max = waveform.adc_max; 
        }
        
        // get waveform charge
        const double charge =  waveform.adc_max - waveform.pedestal;
        assert( charge > 0 );
        
        // update cluster position
        const auto strip_begin = m_geometry_container->get_strip_begin( waveform.layer, waveform.tile, waveform.strip );
        const auto strip_end = m_geometry_container->get_strip_end( waveform.layer, waveform.tile, waveform.strip );
        
        cluster.begin += charge*strip_begin;
        cluster.end += charge*strip_end;
        
        // update cluster charge
        cluster.charge += charge;

      }
      
      cluster.begin *= 1./cluster.charge;
      cluster.end *= 1./cluster.charge;
      cluster.center = 0.5*(cluster.begin+cluster.end);

      cluster.region = cluster.strip/64;
      out.push_back( cluster );
    }
    
    return out;
  }
  
  //_________________________________________________________________
  bool is_signal( const MicromegasRawDataEvaluation::Waveform& waveform, const sample_window_t& sample_window )
  {
    return 
      waveform.rms > 0 &&
      waveform.sample_max >= sample_window.first &&
      waveform.sample_max < sample_window.second &&
      waveform.adc_max >= waveform.pedestal+n_sigma*waveform.rms;
  }
  
  //_________________________________________________________________
  bool is_signal( const MicromegasRawDataEvaluation::Waveform& waveform )
  { return is_signal( waveform, sample_window_signal ); }

  //_________________________________________________________________
  bool is_background( const MicromegasRawDataEvaluation::Waveform& waveform )
  { return is_signal( waveform, sample_window_background ); }
  
  //_________________________________________________________________
  bool accept_event() 
  {
    bool accepted = false;
    
    // loop over tiles
    for( int itile = 0; itile < 8; ++itile )
    {
      const int detid_phi = itile;
      const int detid_z = itile + 8;
      const auto n_clusters_phi = m_fullEventContainer->n_detector_clusters[detid_phi];
      const auto n_clusters_z = m_fullEventContainer->n_detector_clusters[detid_z];
      if( n_clusters_phi >= 1 &&  n_clusters_z >= 1 )
      {
        accepted = true;
        break;
      }
    }    
    return accepted;
  }
  
  //_________________________________________________________________
  void process_event( uint64_t lvl1_bco, const MicromegasRawDataEvaluation::Waveform::List& waveforms )
  {

    // reset full event container
    m_fullEventContainer->Reset();
    m_fullEventContainer->lvl1_bco = lvl1_bco;
    m_fullEventContainer->n_waveforms_all = waveforms.size();
      
    waveform_map_t waveform_map;
    for( const auto& waveform:waveforms )
    {
      if( !waveform.layer ) continue;

//       std::cout << "process_event -"
//         << " layer: " << waveform.layer << " tile: " << waveform.tile << " strip: " << waveform.strip 
//         << " adc_max: " << waveform.adc_max
//         << " pedestal: " << waveform.pedestal
//         << " rms: " << waveform.rms
//         << " signal: " << is_signal( waveform )
//         << std::endl;

      const int detid = waveform.tile + 8*(waveform.layer-55 );
      if( is_signal( waveform ) )
      {        
        waveform_map[detid].insert( waveform );
        ++m_fullEventContainer->n_waveforms_signal;
      }
    }
        
    std::cout << "process_event - bco: 0x" << std::hex << lvl1_bco << std::dec << " signal: " << m_fullEventContainer->n_waveforms_signal << std::endl;
    
    // assign proper size to arrays
    m_fullEventContainer->clusters.clear();
    m_fullEventContainer->n_detector_clusters.resize(16,0);
    m_fullEventContainer->n_region_clusters.resize(64,0);
    m_fullEventContainer->min_cluster_size.resize(16,0);
    m_fullEventContainer->min_cluster_charge.resize(16,0);
    m_fullEventContainer->first_cluster_strip.resize(16,0);

    m_fullEventContainer->n_clusters = 0;
    for( const auto& [detid, waveform_set]:waveform_map )
    {
      auto clusters = get_clusters( waveform_set );
      for( auto&& cluster:clusters ) cluster.is_signal = true;
      
      // process clusters
      m_fullEventContainer->n_detector_clusters[detid] = clusters.size();
      m_fullEventContainer->n_clusters += clusters.size();
      
      unsigned short min_cluster_size = 0;
      double min_cluster_charge = 0;
      
      for( const auto& cluster:clusters )
      { 
        // fill per region counts
        int region_id = cluster.region+4*detid;
        ++m_fullEventContainer->n_region_clusters[region_id];
        
        if( min_cluster_size == 0 || cluster.size < min_cluster_size ) { min_cluster_size = cluster.size; }
        if( min_cluster_charge == 0 || cluster.charge < min_cluster_charge ) { min_cluster_charge = cluster.charge; }                
      }
      
      m_fullEventContainer->min_cluster_size[detid] = min_cluster_size;
      m_fullEventContainer->min_cluster_charge[detid] = min_cluster_charge;
      
      if( !clusters.empty() ) 
      { m_fullEventContainer->first_cluster_strip[detid] = clusters.front().strip; }
      
      // per view clusters
      if( detid < 8 ) m_fullEventContainer->n_phi_clusters += clusters.size();
      else m_fullEventContainer->n_z_clusters += clusters.size();

      // copy into global list
      std::copy( clusters.begin(), clusters.end(), std::back_inserter(m_fullEventContainer->clusters));

    }

    const bool accepted =  !m_apply_event_filter || accept_event();    
    
    if( false )
    {
      std::cout
        << "process_event -"
        << " lvl1_bco: " << lvl1_bco 
        << " waveforms: " << waveforms.size() 
        << " accepted: " << accepted
        << std::endl;
    }

    if( accepted ) 
    {
      m_filtered_lvl1_bco.push_back( lvl1_bco );
      tree_out->Fill(); 
    }
    
    // update counters
    ++m_bco_count;
    if( accepted ) ++m_bco_count_filtered;
    
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

//_____________________________________________________________________________
TString RawDataClusterTree( const TString& inputfilename, const TString& rootfilename )
{
  std::cout << "RawDataClusterTree - inputfilename: " << inputfilename << std::endl;
  std::cout << "RawDataClusterTree - rootfilename: " << rootfilename << std::endl;

  // load geometry
  load_geometry( "micromegas_geometry.root" );

  // filemanager
  FileManager fileManager( inputfilename );
  auto tree = fileManager.GetChain( "T" );
  if( !tree )
  {
    std::cout << "RawDataClusterTree - invalid file" << std::endl;
    return rootfilename;
  }
  
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // output tree
  create_tree( rootfilename );
  
  // keep track of all waveforms associated to a given bco   
  MicromegasRawDataEvaluation::Waveform::List waveforms;
  
  // previous bco and counter
  uint64_t prev_lvl1_bco = 0;
  
  // loop over tree entries
  // const int entries = tree->GetEntries();
  const int entries = 5;
  std::cout << "RawDataClusterTree - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%500) )
    { std::cout << "RawDataClusterTree - entry: " << i << std::endl; }
    
    tree->GetEntry(i);    
    for( const auto& waveform:container->waveforms )
    {
      // skip waveform for which lvl1_bco is zero: they could not be properly associated
      if( !waveform.lvl1_bco ) continue;
      
      // get current bco
      if( std::abs( ((int64_t)waveform.lvl1_bco) - ((int64_t)prev_lvl1_bco) ) > 20 )
      {
        
        // process current waveforms
        if( !waveforms.empty() ) process_event( prev_lvl1_bco, waveforms );
        waveforms.clear();
        
        // update previous bco, store current waveform
        prev_lvl1_bco = waveform.lvl1_bco;
      }
      
      waveforms.push_back( waveform );

    }    
  }

  // process last event 
  if( !waveforms.empty() ) process_event( prev_lvl1_bco, waveforms );

  // save output tree
  save_tree();
  
  if( m_apply_event_filter )
  {
    print_list( "lvl1_bco", "uint64_t", m_filtered_lvl1_bco );    
    std::cout << "RawDataClusterTree -"
      << " m_bco_count: " << m_bco_count
      << " m_bco_count_filtered: " << m_bco_count_filtered
      << " ratio: " << (double)m_bco_count_filtered/m_bco_count 
      << std::endl;
  }
    
  std::cout << "done." << std::endl;
  
  return rootfilename;
}

//_____________________________________________________________________________
TString RawDataClusterTree( int runNumber = 20121 )
{
  
  const TString inputfilename = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-0000-test.root", runNumber );
  // const TString inputfilename = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString rootfilename = m_apply_event_filter ? 
    Form( "Rootfiles/RawDataClusterTree-%08i-0000-test_filtered.root", runNumber ):
    Form( "Rootfiles/RawDataClusterTree-%08i-0000-test.root", runNumber );
 
  return RawDataClusterTree( inputfilename, rootfilename );
  
}
