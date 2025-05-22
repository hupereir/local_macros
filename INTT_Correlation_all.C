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
  class InttEvent : public TObject
  {
    public:
    int evtSeq = 0;
    Long64_t bco = 0;
    int fNhits = 0;
    TClonesArray* fhitArray = nullptr;
    ClassDef(InttEvent, 2)

  };

  static constexpr int n_intt_boards = 8;
  std::array<InttEvent*,n_intt_boards> intt_event_array;

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
  FullEventContainer* tpot_container = nullptr;
  
  //_________________________________________________
  using offset_t = std::pair<size_t, int>;
  std::optional<offset_t> get_intt_event_offset( TTree* intt_tree, InttEvent* intt_event, TTree* tpot_tree )
  {

    // store first 100 deltas from tpot tree
    const int max_entry = 100;
    int tpot_entries = std::min<int>( max_entry, tpot_tree->GetEntries() );

    std::vector<uint64_t> tpot_bco;
    for( int i = 0; i < tpot_entries; ++i )
    {
      tpot_tree->GetEntry(i);
      tpot_bco.push_back( tpot_container->lvl1_bco );
    }

    // store first 100 deltas from intt tree
    int intt_entries = std::min<int>( max_entry, intt_tree->GetEntries() );

    std::vector<uint64_t> intt_bco;
    const int intt_base_offset = 0;
    // const int intt_base_offset = 560;
    for( int i = 0; i < intt_entries; ++i )
    {
      intt_tree->GetEntry(i+intt_base_offset);
      intt_bco.push_back( intt_event->bco );
    }
     
    // print
    if( false )
    {
      for( size_t i = 0; i < std::min( tpot_bco.size(), intt_bco.size() ); ++i )
      { std::cout << "get_intt_event_offset - entry: " << i << " tpot_bco: " << tpot_bco[i] << " intt_bco: " << intt_bco[i] << std::endl; }
    }
       
    bool found = false;
    int first_tpot_entry = 0;
    int intt_offset = 0;
    for( size_t i = 0; i < intt_bco.size()-1; ++i )
    {
      for( size_t itpot = 0; itpot < tpot_bco.size() && !found; ++itpot )
      {
        if( std::abs<int64_t>( int64_t(tpot_bco[itpot])-intt_bco[i] ) <= 1 )
        {
          first_tpot_entry = itpot;
          intt_offset = int(i)-itpot;
          found = true;
        }
      }      
    }
    
    if( !found ) std::cout << "get_intt_event_offset - no matching events found" << std::endl;
    else std::cout 
      << "get_intt_event_offset -"
      << " first_tpot_entry: " << first_tpot_entry
      << " intt_offset: " << intt_offset+intt_base_offset
      << std::endl;
    
    return found ? std::optional(std::make_pair( first_tpot_entry, intt_offset+intt_base_offset )):std::nullopt;
  }
  
}

//_________________________________________________
void INTT_Correlation_all( int runnumber = 20446 )
{
  set_style( false );

  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  const TString pdffilename = Form( "Figures/INTT_Correlation-%08i-0000-all.pdf", runnumber );
  const TString pngfilename = Form( "Figures/INTT_Correlation-%08i-0000-all.png", runnumber );
  const TString rootfilename = Form( "Rootfiles/INTT_Correlation-%08i-0000-all.root", runnumber );

  std::cout << "INTT_Correlation - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  PdfDocument pdfDocument( pdffilename );
  RootFile rootfile( rootfilename );

  // read all INTT input files, setup tree
  std::array<TTree*, n_intt_boards> intt_tree_array;
  for( int i = 0 ; i < n_intt_boards; ++i )
  {
    const TString intt_inputfilename = Form( "/phenix/u/hpereira/sphenix/work/intt/LUSTRE/beam_intt%i-%08i-0000_event_base.root", i, runnumber );
    std::cout << "INTT_Correlation - intt_inputfilename: " << intt_inputfilename << std::endl;
    
    auto intt_tfile = TFile::Open( intt_inputfilename, "READ" );
    intt_tree_array[i] = static_cast<TTree*>( intt_tfile->Get("tree") );
    intt_event_array[i] = nullptr;
    intt_tree_array[i]->SetBranchAddress( "event", &intt_event_array[i] );
  }
  
  // read TPOT input file, setup tree
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  tpot_tree->SetBranchAddress( "Event", &tpot_container );

  // find event offsets
  size_t first_tpot_entry = 0;
  std::array<int, n_intt_boards> intt_offset_array;
  for( int i = 0 ; i < n_intt_boards; ++i )
  {  
    const auto& intt_tree = intt_tree_array[i];
    const auto& intt_event = intt_event_array[i];
    auto offset = get_intt_event_offset( intt_tree, intt_event, tpot_tree );
    if( !offset ) return;
    first_tpot_entry = std::max( first_tpot_entry, offset->first );
    intt_offset_array[i] = offset->second;
  }
    
  // define first entry, including offset
  size_t tpot_entry = first_tpot_entry;
  
  std::array<size_t, n_intt_boards> intt_entry_array;
  for( int i = 0; i < n_intt_boards; ++i )
  { intt_entry_array[i] = first_tpot_entry+intt_offset_array[i]; }

  const int max_entry = 1e6;
  int tpot_entries = std::min<int>( max_entry, tpot_tree->GetEntries() );

  std::array<int,n_intt_boards> intt_entries_array;
  for( int i = 0; i < n_intt_boards; ++i )
  { intt_entries_array[i] = std::min<int>( max_entry, intt_tree_array[i]->GetEntries() ); }

  // correlation histogram
  auto h_correlation = new TH2F( "h_correlation", "", 320, 0, 320, 320, 0, 2200*n_intt_boards );
  h_correlation->GetXaxis()->SetTitle( "TPOT Clusters" );
  h_correlation->GetYaxis()->SetTitle( "INTT Hits" );
  h_correlation->GetYaxis()->SetTitleOffset( 1.4 );
  rootfile.Add( h_correlation );
  
  bool empty = false;
  while( !empty )
  {

    tpot_tree->GetEntry( tpot_entry );
    auto tpot_bco = tpot_container->lvl1_bco;
    
    // we assume that all INTT are synchronized
    intt_tree_array[0]->GetEntry( intt_entry_array[0] );
    auto intt_bco = intt_event_array[0]->bco;

    std::cout << "TPOT_Correlation -"
      << " tpot_entry: " << tpot_entry
      << " intt_entry: " << intt_entry_array[0]
      << " tpot bco: " << tpot_bco
      << " intt bco: " << intt_bco
      << std::endl;

    /* 
     * compare intt time and tpot time. 
     * The latter can be significantly bigger than the former when for some reason tpot events are dropper 
     * one needs to add some fuzziness to the comparison. 1000 seems a good number
     */
    bool corrected = false;
    while( tpot_bco > intt_bco + 1 && intt_entry_array[0]<intt_entries_array[0] )
    {      
      corrected = true;
      
      for( int i = 0; i < n_intt_boards; ++i )
      { ++intt_entry_array[i]; }
      
      
      intt_tree_array[0]->GetEntry( intt_entry_array[0] );
      intt_bco = intt_event_array[0]->bco;
      
      // printout
      std::cout << "TPOT_Correlation -"
        << " tpot_entry: " << tpot_entry
        << " intt_entry: " << intt_entry_array[0]
        << " tpot bco: " << tpot_bco
        << " intt bco: " << intt_bco
        << " - corrected (TPOT)"
        << std::endl;
    }

    while( intt_bco > tpot_bco + 1 && tpot_entry<tpot_entries )
    {
      
      corrected = true;
      ++tpot_entry;
      tpot_tree->GetEntry( tpot_entry );
      tpot_bco = tpot_container->lvl1_bco;
      
      // printout
      std::cout << "TPOT_Correlation -"
        << " tpot_entry: " << tpot_entry
        << " intt_entry: " << intt_entry_array[0]
        << " tpot bco: " << tpot_bco
        << " intt bco: " << intt_bco
        << " - corrected (INTT)"
        << std::endl;
    }
    
    if( !corrected ) 
    {
      int nhits = 0;
      
      // compute total number of hits
      for( int i = 0; i < n_intt_boards; ++i )
      {
        intt_tree_array[i]->GetEntry( intt_entry_array[i] );
        nhits += intt_event_array[i]->fNhits;
        assert(  intt_event_array[i] ==  intt_event_array[0] );
      }
      
      // fill correlation histogram
      h_correlation->Fill( tpot_container->n_clusters, nhits );
      
    }

    // increment entries
    if( ++tpot_entry >= tpot_entries ) empty = true;
    for( int i = 0; i < n_intt_boards; ++i )
    {  if( ++intt_entry_array[i] >= intt_entries_array[i] ) empty = true; }
    
  }

  if( true )
  {
    TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
    h_correlation->GetXaxis()->SetTitle( "TPOT clusters" );
    h_correlation->GetYaxis()->SetTitle( "INTT hits (all)" );
    h_correlation->Draw("colz");
    h_correlation->GetYaxis()->SetTitleOffset( 1.6 );
    gPad->SetLeftMargin( 0.16);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
      
    auto text = new TPaveText(0.57,0.16,0.86,0.31, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(1010);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText("#it{#bf{sPHENIX}} Internal");
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText("07/11/2023");
    text->Draw();
    
    pdfDocument.Add(cv);
    cv->SaveAs( pngfilename );
  } 
    
}

void process_all()
{
  for( int runnumber:{ 20445, 20446 } )
  { INTT_Correlation_all(runnumber); }
}
