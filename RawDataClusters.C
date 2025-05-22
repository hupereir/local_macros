#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>

#include <memory>
#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

#include <boost/container_hash/hash.hpp>

bool operator < (const MicromegasRawDataEvaluation::Waveform& lhs, const MicromegasRawDataEvaluation::Waveform& rhs )
{ return lhs.strip < rhs.strip; }

namespace
{
  
  class cluster_t: public std::vector<MicromegasRawDataEvaluation::Waveform>
  {
    
    public:
    bool is_signal = true;
    
  };
  
  // sort waveforms per detid
  using cluster_list_t = std::vector<cluster_t>;
  using waveform_set_t = std::set<MicromegasRawDataEvaluation::Waveform>;
  using waveform_map_t = std::map<int,waveform_set_t>;
  
  // number of sigma above pedestal to remove noise
  static constexpr double n_sigma = 5.;
  
  // detector names
  std::vector<std::string> detector_names;
  
  // charge vs size vs region id
  TH3* h_cluster_charge = nullptr;
  
  // charge vs size vs region id
  TH3* h_cluster_charge_background = nullptr;
  
  // multiplicity vs region
  TH2* h_cluster_mult = nullptr;

  // phi/z correlation
  TH3* h_cluster_mult_correlation = nullptr;
  
  // phi/z correlation
  TH2* h_cluster_mult_correlation_all = nullptr;

  MicromegasMapping mapping;

  MicromegasCalibrationData calibration_data;
   
  // sample windows to determine signal and background clusters
  using sample_window_t = std::pair<unsigned short, unsigned short>;
  // sample_window_t sample_window_signal = {20, 45};
  sample_window_t sample_window_signal = {15, 40};
  sample_window_t sample_window_background = {5, 15 };
  
  // save detector names
  void save_detector_names()
  {
    
    // get detector names that match tile/layer
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name ) );
    }
  }
  
  // create histograms
  void create_histograms()
  {
    static constexpr double charge_max = 4000;
    static constexpr int n_charge_bins = 500;
    h_cluster_charge = new TH3I( "h_cluster_charge", "h_cluster_charge", 64, 0, 64, 20, 0, 20, n_charge_bins, 0, charge_max );
    h_cluster_charge->GetXaxis()->SetTitle( "region id" );
    h_cluster_charge->GetYaxis()->SetTitle( "cluster size" );
    h_cluster_charge->GetZaxis()->SetTitle( "cluster charge (adc)" );
    
    h_cluster_charge_background = new TH3I( "h_cluster_charge_background", "h_cluster_charge_background", 64, 0, 64, 20, 0, 20, n_charge_bins, 0, charge_max );
    h_cluster_charge_background->GetXaxis()->SetTitle( "region id" );
    h_cluster_charge_background->GetYaxis()->SetTitle( "cluster size" );
    h_cluster_charge_background->GetZaxis()->SetTitle( "cluster charge (adc)" );

    h_cluster_mult = new TH2I( "h_cluster_mult", "h_cluster_mult", 16, 0, 16, 100, 0, 100 );
    h_cluster_mult->GetXaxis()->SetTitle( "detector id" );
    h_cluster_mult->GetYaxis()->SetTitle( "cluster multiplicity" );

    h_cluster_mult_correlation = new TH3I( "h_cluster_mult_correlation", "h_cluster_mult_correlation", 8, 0, 8, 100, 0, 100, 100, 0, 100 );
    h_cluster_mult_correlation->GetXaxis()->SetTitle( "detector id" );
    h_cluster_mult_correlation->GetYaxis()->SetTitle( "cluster multiplicity (phi)" );
    h_cluster_mult_correlation->GetZaxis()->SetTitle( "cluster multiplicity (z)" );

    h_cluster_mult_correlation_all = new TH2I( "h_cluster_mult_correlation_all", "h_cluster_mult_correlation_all", 100, 0, 100, 100, 0, 100 );
    h_cluster_mult_correlation_all->GetXaxis()->SetTitle( "cluster multiplicity (phi)" );
    h_cluster_mult_correlation_all->GetYaxis()->SetTitle( "cluster multiplicity (z)" );

  }
    
  // cluster processing
  //_________________________________________________________________
  void process_cluster( int detid, const cluster_t& cluster )
  {
    double charge_max = 0;
    unsigned short strip_max = 0;
    double charge_total = 0;
    for( const auto& waveform:cluster )
    { 
      const double charge = waveform.adc_max-waveform.pedestal;
      if( charge > charge_max ) 
      {
        charge_max = charge;
        strip_max = waveform.strip;
      }
      charge_total += charge;
    }
    
    // get region
    unsigned short iregion = strip_max / 64;
    int region_id = iregion + 4*detid;

    // fill
    if( cluster.is_signal ) h_cluster_charge->Fill(region_id,cluster.size(),charge_total);
    else h_cluster_charge_background->Fill(region_id,cluster.size(),charge_total); 
  }

  //_________________________________________________________________
  void process_clusters( int detid, const cluster_list_t clusters )
  { 
    for( const auto& cluster:clusters ) 
    { process_cluster( detid, cluster ); }  
  }
  
  //_________________________________________________________________
  cluster_list_t get_clusters( const waveform_set_t& waveform_set )
  {
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
    
    return clusters;
  }

  //_________________________________________________________________
  bool is_signal( const MicromegasRawDataEvaluation::Waveform& waveform )
  {
    return 
      waveform.rms > 0 &&
      waveform.sample_max >= sample_window_signal.first &&
      waveform.sample_max < sample_window_signal.second &&
      waveform.adc_max > waveform.pedestal+n_sigma*waveform.rms;
  }

  //_________________________________________________________________
  bool is_background( const MicromegasRawDataEvaluation::Waveform& waveform )
  {
    return 
      waveform.rms > 0 &&
      waveform.sample_max >= sample_window_background.first &&
      waveform.sample_max < sample_window_background.second &&
      waveform.adc_max > waveform.pedestal+n_sigma*waveform.rms;
  }
  
  //_________________________________________________________________
  void process_event( uint64_t /*lvl1_bco*/, const MicromegasRawDataEvaluation::Waveform::List& waveforms )
  {

    waveform_map_t waveform_map;
    waveform_map_t waveform_map_background;

    for( const auto& waveform:waveforms )
    {
      if( !waveform.layer ) continue;
      const int detid = waveform.tile + 8*(waveform.layer-55 );
      if( is_signal( waveform ) ) waveform_map[detid].insert( waveform );
      if( is_background( waveform ) ) waveform_map_background[detid].insert( waveform );
    }
        
    // store cluster multiplicity
    std::array<size_t,16> cluster_mult = {};
    for( const auto& [detid, waveform_set]:waveform_map )
    {
      auto clusters = get_clusters( waveform_set );
      for( auto&& cluster:clusters ) cluster.is_signal = true;
      
      // process clusters
      cluster_mult[detid] = clusters.size();
      process_clusters( detid, clusters );
    }    
    
    // fill cluster multiplicity
    for(int detid=0; detid<16;++detid)
    { h_cluster_mult->Fill(detid, cluster_mult[detid]); }
    
    // fill cluster correlations
    int mult_phi_all = 0;
    int mult_z_all = 0;
    for( int tile =0; tile < 8; ++tile )
    {
      const auto mult_phi = cluster_mult[tile];
      const auto mult_z = cluster_mult[tile+8];
      
      
      // std::cout << "process_event - correlation: " << detector_names[tile] << " " << mult_phi << std::endl;
      // std::cout << "process_event - correlation: " << detector_names[tile+8] << " " << mult_z << std::endl;
      if( tile > 0 )
      {
        // ignore tile 0 because of disconnected FEE
        mult_phi_all += mult_phi;
        mult_z_all += mult_z;
      }
      h_cluster_mult_correlation->Fill( tile, mult_phi, mult_z );
    }
    
    h_cluster_mult_correlation_all->Fill( mult_phi_all, mult_z_all );
    
    // process background clusters
    for( const auto& [detid, waveform_set]:waveform_map_background )
    {
      auto clusters = get_clusters( waveform_set );
      for( auto&& cluster:clusters ) cluster.is_signal = false;
      
      // process clusters
      process_clusters( detid, clusters );
    }    
  }
  
  //_________________________________________________________________
  void draw_histograms( PdfDocument& pdfDocument )
  {

    if( true )
    {
      for( int ilayer = 0; ilayer <2; ++ilayer )
        for( int itile = 0; itile <8; ++itile )
      {
        const int detid = itile + 8*ilayer;
        auto cvname = Form( "cv_%i", detid);
        auto cv = std::make_unique<TCanvas>( cvname, cvname, 900,500 );
        cv->Divide( 2, 2 );
        for( int iregion = 0; iregion<4; ++iregion )
        {
          cv->cd( iregion+1);
          int region_id = iregion +4*detid;
          {
            h_cluster_charge->GetXaxis()->SetRange(region_id+1, region_id+1);
            auto h = static_cast<TH1*>( h_cluster_charge->Project3D( "z" ) );          
            const auto hname = Form( "h_%s_R%i", detector_names[detid].c_str(), 4-iregion );
            const auto htitle = Form( "charge %s_R%i", detector_names[detid].c_str(), 4-iregion );
            h->SetName( hname );
            h->SetTitle( htitle );
            h->SetLineColor( 2 );
            h->Draw();
          }
          
          gPad->SetLogy(true);
          
          // background histogram
          {
            h_cluster_charge_background->GetXaxis()->SetRange(region_id+1, region_id+1);
            auto h = static_cast<TH1*>( h_cluster_charge_background->Project3D( "z" ) );          
            const auto hname = Form( "h_%s_R%i_bg", detector_names[detid].c_str(), 4-iregion );
            const auto htitle = Form( "charge %s_R%i (bg)", detector_names[detid].c_str(), 4-iregion );
            h->SetName( hname );
            h->SetTitle( htitle );
            h->Draw( "same" );
          }
        }
        pdfDocument.Add( cv.get() );
      }
    }
  }
  
  //_________________________________________________________________
  void save_histograms( RootFile& rootfile )
  {
    rootfile.Add( h_cluster_charge );
    rootfile.Add( h_cluster_charge_background );
    rootfile.Add( h_cluster_mult );
    rootfile.Add( h_cluster_mult_correlation );
    rootfile.Add( h_cluster_mult_correlation_all );
  }
}

//_____________________________________________________________________________
TString RawDataClusters( int runNumber = 21161 )
{
  const TString inputfilename = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-0000-full.root", runNumber );
  const TString pdffilename = Form( "Figures/RawDataClusters-%08i-0000.pdf", runNumber );
  const TString rootfilename = Form( "Rootfiles/RawDataClusters-%08i-0000.root", runNumber );
  
  std::cout << "RawDataClusters - inputfilename: " << inputfilename << std::endl;
  std::cout << "RawDataClusters - pdffilename: " << pdffilename << std::endl;
  std::cout << "RawDataClusters - rootfilename: " << rootfilename << std::endl;

  save_detector_names();
  create_histograms();
  
  // ouput pdffile
  PdfDocument pdfDocument( pdffilename );

  // output rootfile
  RootFile rootFile( rootfilename );

  // filemanager
  FileManager fileManager( inputfilename );
  auto tree = fileManager.GetChain( "T" );
  if( !tree )
  {
    std::cout << "RawDataClusters - invalid file" << std::endl;
    return pdffilename;
  }
  
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // keep track of all waveforms associated to a given bco 
  using waveform_map_t = std::map<uint64_t, MicromegasRawDataEvaluation::Waveform::List>;
  waveform_map_t waveform_map;
  
  // loop over tree entries
  const int entries = tree->GetEntries();
  std::cout << "RawDataClusters - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    if( !(i%500) ) 
    { std::cout << "RawDataClusters - entry: " << i << std::endl; }
    
    tree->GetEntry(i);
        
    if( false )
    {
      std::cout << "RawDataClusters -"
        << " lvl1_bco_list: " << container->lvl1_bco_list.size()
        << " " <<  container->lvl1_bco_list[0]
        << " " <<  container->lvl1_bco_list[1]
        << std::endl;
      
      std::cout << "RawDataClusters -"
        << " lvl1_count_list: " << container->lvl1_count_list.size()
        << " " <<  container->lvl1_count_list[0]
        << " " <<  container->lvl1_count_list[1]
        << std::endl;
    }
    
    // loop over waveforms
    for( const auto& waveform:container->waveforms )
    { waveform_map[waveform.lvl1_bco].push_back( waveform ); }    
  }

  for( const auto& [lvl1_bco,waveforms]:waveform_map)
  {
    std::cout
      << "RawDataClusters -"
      << " bco: " << lvl1_bco 
      << " waveforms: " << waveforms.size() 
      << std::endl;
        
    // make plot
    process_event( lvl1_bco, waveforms );
  }

  // save histograms
  draw_histograms( pdfDocument );
  save_histograms( rootFile );
  std::cout << "done." << std::endl;
  
  return pdffilename;
}

//__________________________________________________________________
void process_all()
{
  using run_map_t = std::vector<std::pair<int,int>>;
//   run_map_t run_map = 
//   {
//     { 50, 20464 },
//     { 100, 20465 },
//     { 150, 20466 },
//     { 200, 20467 },
//     { 250, 20468 },
//     { 300, 20469 },
//     { 350, 20470 },
//     { 400, 20471 }
//   };
  
  run_map_t run_map = 
  {
//     { 430, 21079 },
//     { 420, 21080 },
//     { 410, 21081 },
//     { 400, 21082 },
//     { 390, 21083 },
//     { 380, 21084 },
//     { 370, 21153 },
//     { 360, 21154 },
//     { 350, 21155 },
//     { 340, 21156 },
//     { 330, 21157 },
//     { 320, 21158 },
//     
//     { 430, 21159 },
//     { 420, 21160 },
    { 410, 21161 },
    { 400, 21162 },
    { 390, 21163 },
    { 380, 21164 },
    { 370, 21166 },
    { 360, 21168 },
    { 350, 21169 },
    { 340, 21170 },
    { 330, 21171 },
    { 320, 21172 }
  };
  
  for( const auto& [hv, runnumber]:run_map )
  {  RawDataClusters(runnumber);  }
}
