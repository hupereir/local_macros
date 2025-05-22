#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include <optional>

namespace 
{
  
  //____________________________________________________________________________
  class MBDContainer
  {
    public:
    int   evt = 0;
    unsigned short clk = 0;
    unsigned short femclk = 0;
    float bqs = 0;
    float bqn = 0;
  };
  
  auto mbd_container = std::make_unique<MBDContainer>();
  static constexpr uint64_t mbd_max_clk = 0x10000;

  //_______________________________________________
  void setup_mbd_tree( TTree* mbd_tree )
  {    
    mbd_tree->SetBranchAddress("evt",&mbd_container.get()->evt);
    mbd_tree->SetBranchAddress("clk",&mbd_container.get()->clk);
    mbd_tree->SetBranchAddress("femclk",&mbd_container.get()->femclk);
    mbd_tree->SetBranchAddress("bqs",&mbd_container.get()->bqs);
    mbd_tree->SetBranchAddress("bqn",&mbd_container.get()->bqn);
  }
    
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
      n_detector_clusters.clear();
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
    
    // number of clusters per detector
    std::vector<unsigned short> n_detector_clusters;
    
    // number of clusters
    unsigned short n_phi_clusters = 0;
    
    // number of clusters
    unsigned short n_z_clusters = 0;
    
    ClassDef(FullEventContainer,1)
      
  };
  
  // auto tpot_container = std::make_unique<FullEventContainer>();
  auto tpot_container = new FullEventContainer;

  //_______________________________________________
  void setup_tpot_tree( TTree* tpot_tree )
  {    
    tpot_tree->SetBranchAddress( "Event", &tpot_container );
  }

  //_______________________________________________
  uint64_t get_mbd_delta( uint64_t first, uint64_t second )
  {
    while( second < first ) second += mbd_max_clk;
    return second - first;
  }
  
  //_________________________________________________
  using offset_t = std::pair<size_t, int>;
  std::optional<offset_t> get_mbd_event_offset( TTree* mbd_tree, TTree* tpot_tree )
  {

    // store first 100 deltas from tpot tree
    const int max_entry = 100;
    int tpot_entries = std::min<int>( max_entry, tpot_tree->GetEntries() );

    std::vector<uint64_t> tpot_deltas;
    uint64_t tpot_prev_bco = 0;
    for( int i = 0; i < tpot_entries; ++i )
    {
      tpot_tree->GetEntry(i);
      tpot_deltas.push_back( tpot_container->lvl1_bco - tpot_prev_bco );
      tpot_prev_bco = tpot_container->lvl1_bco;
    }

    // store first 100 deltas from mbd tree
    int mbd_entries = std::min<int>( max_entry, mbd_tree->GetEntries() );

    std::vector<uint64_t> mbd_deltas;
    uint64_t mbd_prev_clk = 0;
    for( int i = 0; i < mbd_entries; ++i )
    {
      mbd_tree->GetEntry(i);
      mbd_deltas.push_back( get_mbd_delta( mbd_prev_clk, mbd_container->clk ));
      mbd_prev_clk = mbd_container->clk;
    }
    
    // print
    if( false )
    {
      for( size_t i = 0; i < std::min( tpot_deltas.size(), mbd_deltas.size() ); ++i )
      { std::cout << "get_mbd_event_offset - entry: " << i << " tpot_delta: " << tpot_deltas[i] << " mbd_delta: " << mbd_deltas[i] << std::endl; }
    }
    
    bool found = false;
    int first_tpot_entry = 0;
    int mbd_offset = 0;
    for( size_t i = 0; i < mbd_deltas.size()-1; ++i )
    {
      for( size_t itpot = 0; itpot < tpot_deltas.size()-1 && !found; ++itpot )
      {
        if( tpot_deltas[itpot] == mbd_deltas[i] && tpot_deltas[itpot+1] == mbd_deltas[i+1] )
        {
          first_tpot_entry = itpot;
          mbd_offset = int(i)-itpot;
          found = true;
        }
      }      
    }
    
    if( !found ) std::cout << "get_mbd_event_offset - no matching events found" << std::endl;
    else std::cout 
      << "get_mbd_event_offset -"
      << " first_tpot_entry: " << first_tpot_entry
      << " mbd_offset: " << mbd_offset
      << std::endl;
    
    return found ? std::optional(std::make_pair( first_tpot_entry, mbd_offset )):std::nullopt;
  }
  
}

//_________________________________________________
void MBD_Synchronization()
{
  
  static constexpr int NUM_ARMS = 2;

  const int runnumber = 20445;
  const TString mbd_inputfilename = Form( "/sphenix/user/chiu/sphenix_bbc/run2023/beam_seb18-%08i-0000_mbd.root", runnumber );
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000.root", runnumber );
  
  std::cout << "MBD_Synchronization - mbd_inputfilename: " << mbd_inputfilename << std::endl;
  std::cout << "MBD_Synchronization - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  
  auto mbd_tfile = TFile::Open( mbd_inputfilename, "READ" );
  auto mbd_tree = static_cast<TTree*>( mbd_tfile->Get("t") );
  setup_mbd_tree( mbd_tree );
  
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  setup_tpot_tree( tpot_tree );
  
  auto offset = get_mbd_event_offset( mbd_tree, tpot_tree );
  if( !offset ) return;
  const auto& [first_tpot_entry, mbd_offset] = *offset;
  
  // loop over entries using offset
  size_t tpot_entry = first_tpot_entry;
  size_t mbd_entry = first_tpot_entry+mbd_offset;
  
  uint64_t tpot_prev_bco = 0;
  uint64_t mbd_prev_clk = 0;
  
  // store first 100 deltas from tpot tree
  const int max_entry = 10000;
  int tpot_entries = std::min<int>( max_entry, tpot_tree->GetEntries() );
  int mbd_entries = std::min<int>( max_entry, mbd_tree->GetEntries() );

  bool first = true;
  while( tpot_entry < tpot_entries && mbd_entry < mbd_entries )
  {
    tpot_tree->GetEntry( tpot_entry );
    uint64_t tpot_delta = tpot_container->lvl1_bco-tpot_prev_bco;
    tpot_prev_bco = tpot_container->lvl1_bco;

    mbd_tree->GetEntry( mbd_entry );
    auto mbd_delta = get_mbd_delta( mbd_prev_clk, mbd_container->clk );
    mbd_prev_clk = mbd_container->clk;

    std::cout << "MBD_Synchronization -"
      << " tpot_entry: " << tpot_entry
      << " mbd_entry: " << mbd_entry 
      << " tpot_delta: " << tpot_delta
      << " mbd_delta: " << mbd_delta
      << std::endl;
    
    while( !first && (tpot_delta  > mbd_delta + 1000) )
    {
      ++mbd_entry;
      mbd_tree->GetEntry( mbd_entry );
      mbd_delta += get_mbd_delta( mbd_prev_clk, mbd_container->clk );
      mbd_prev_clk = mbd_container->clk;
      std::cout << "MBD_Synchronization -"
        << " tpot_entry: " << tpot_entry
        << " mbd_entry: " << mbd_entry 
        << " tpot_delta: " << tpot_delta
        << " mbd_delta: " << mbd_delta
        << " correcting" 
        << std::endl;
    }
    
    first = false;
    ++tpot_entry;
    ++mbd_entry;
    
  }
  
  
}
