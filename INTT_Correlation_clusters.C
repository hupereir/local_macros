#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>

#include <optional>
#include <cmath>
#include <cstdint>

/// INTT event structure
class InttEvent : public TObject
{
  public:
  int evtSeq = 0;
  Long64_t bco = 0;
  int fNhits = 0;
  TClonesArray* fhitArray = nullptr;
  int fNclusters = 0;
  ClassDef(InttEvent, 2)
    
};

/// TPOT event structure
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

namespace
{
  
  // running intt event
  InttEvent* intt_event = nullptr;

  // running tpot container
  FullEventContainer* tpot_container = nullptr;

  using offset_t = std::pair<size_t, int>;
  std::optional<offset_t> get_intt_event_offset( TTree* intt_tree, TTree* tpot_tree )
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
    if( true )
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
        if( fabs( int64_t(tpot_bco[itpot])-intt_bco[i] ) <= 1 )
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
void INTT_Correlation_clusters(const int runnumber = 20446)
{
  set_style( false );

  const int intt_id = 0;
  
  const TString intt_inputfilename = Form( "Rootfiles/intt_clusters/beam_intt%i-%08i-0000_clusters.root", intt_id, runnumber );
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  const TString pdffilename = Form( "Figures/INTT_Correlation_clusters-%08i-0000_intt%i.pdf", runnumber, intt_id );
  const TString pngfilename = Form( "Figures/INTT_Correlation_clusters-%08i-0000_intt%i.png", runnumber, intt_id );
  const TString rootfilename = Form( "Rootfiles/INTT_Correlation_clusters-%08i-0000_intt%i.root", runnumber, intt_id );

  std::cout << "INTT_Correlation - intt_inputfilename: " << intt_inputfilename << std::endl;
  std::cout << "INTT_Correlation - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  PdfDocument pdfDocument( pdffilename );
  RootFile rootfile( rootfilename );

  auto intt_tfile = TFile::Open( intt_inputfilename, "READ" );
  auto intt_tree = static_cast<TTree*>( intt_tfile->Get("tree") );
  intt_tree->SetBranchAddress( "event", &intt_event );

  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  tpot_tree->SetBranchAddress( "Event", &tpot_container );

  auto offset = get_intt_event_offset( intt_tree, tpot_tree );
  if( !offset ) return;
  const auto& [first_tpot_entry, intt_offset] = *offset;
    
  // loop over entries using offset
  size_t tpot_entry = first_tpot_entry;
  size_t intt_entry = first_tpot_entry+intt_offset;

  const size_t max_entry = 1e6;
  auto intt_entries = std::min<size_t>( max_entry, intt_tree->GetEntries() );
  auto tpot_entries = std::min<size_t>( max_entry, tpot_tree->GetEntries() );

  // correlation histogram
  auto h_correlation = new TH2F( "h_correlation", "", 320, 0, 320, 320, 0, 1200 );
  h_correlation->GetXaxis()->SetTitle( "TPOT clusters" );
  h_correlation->GetYaxis()->SetTitle( Form( "INTT clusters (INTT%i)", intt_id ) );
  h_correlation->GetYaxis()->SetTitleOffset( 1.6 );
  h_correlation->GetXaxis()->SetNdivisions(505);
  
  // add rootfile
  rootfile.Add( h_correlation );
  
  while( tpot_entry < tpot_entries && intt_entry < intt_entries )
  {
    tpot_tree->GetEntry( tpot_entry );
    auto tpot_bco = tpot_container->lvl1_bco;
    
    intt_tree->GetEntry( intt_entry );
    auto intt_bco = intt_event->bco;

    std::cout << "TPOT_Correlation -"
      << " tpot_entry: " << tpot_entry
      << " intt_entry: " << intt_entry
      << " tpot bco: " << tpot_bco
      << " intt bco: " << intt_bco
      << std::endl;
    
    /* 
     * compare intt time and tpot time. 
     * The latter can be significantly bigger than the former when for some reason tpot events are dropper 
     * one needs to add some fuzziness to the comparison. 1000 seems a good number
     */
    bool corrected = false;
    while( tpot_bco > intt_bco + 1 && intt_entry<intt_entries )
    {
      
      corrected = true;
      ++intt_entry;
      intt_tree->GetEntry( intt_entry );
      intt_bco = intt_event->bco;
      
      // printout
      std::cout << "TPOT_Correlation -"
        << " tpot_entry: " << tpot_entry
        << " intt_entry: " << intt_entry
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
        << " intt_entry: " << intt_entry
        << " tpot bco: " << tpot_bco
        << " intt bco: " << intt_bco
        << " - corrected (INTT)"
        << std::endl;
    }
    
    if( !corrected ) 
    { 
      h_correlation->Fill( tpot_container->n_clusters, intt_event->fNclusters );
    }

    ++tpot_entry;
    ++intt_entry;  
  }
 
  if( true )
  {
    TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
    h_correlation->Draw("colz");
    gPad->SetLeftMargin( 0.16);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
      
    auto text = new TPaveText(0.57,0.16,0.86,0.36, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(1010);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( "" );
    text->AddText( "#it{#bf{sPHENIX}} Internal");
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText("Part of INTT used (INTT0)");
    text->AddText(Form("07/11/2023, Run %i", runnumber ));
    text->Draw();
    
    pdfDocument.Add(cv);
    cv->SaveAs(pngfilename);
  } 
}

//_______________________________________________________-
void process_all()
{
  for( const int runnumber: {20445, 20446} )
  { INTT_Correlation_clusters( runnumber ); }
}
