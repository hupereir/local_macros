#include <RootUtil/PdfDocument.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

void RawDataTestBCO( int runnumber = 24248 )
{
  const TString inputFile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-000?-test.root", runnumber );
  std::cout << "RawDataTestBCO - inputFile: " << inputFile << std::endl;
  
  // input
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // keep track of number of waveform per fee and per bco
  using bco_map_t = std::map<unsigned int, unsigned short>;
  using fee_bco_map_t = std::map<unsigned short, bco_map_t>;
  fee_bco_map_t fee_bco_map;
  
  // loop over tree entries
  const int entries = tree->GetEntries();
  // const int entries = 10000;
  std::cout << "RawDataTestBCO - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%500) )
    { std::cout << "RawDataTestBCO - entry: " << i << std::endl; }
    
    tree->GetEntry(i);    

    // loop over waveforms
    for( const auto& waveform:container->waveforms )
    { ++fee_bco_map[waveform.fee_id][waveform.fee_bco]; }
  }
  
  // make histograms
  using histogram_map_t = std::map<unsigned short, TH1*>;
  histogram_map_t histograms;
  for( const auto& [fee, bco_map]:fee_bco_map )
  {
    auto h = new TH1I( Form( "h_%hu", fee ), Form( "n_waveforms (fee %hu)", fee ), 1024, 0, 1024 );
    histograms[fee] = h;
    for( const auto& [bco, entries]:bco_map )
    { h->Fill( entries ); }
  }

  // draw
  auto cv = new TCanvas( "cv", "", 900, 900 );
  cv->Divide( 4, 4 );
  int cvid = 0;
  for( const auto& [fee, h]:histograms )
  {
    cv->cd( ++cvid );
    h->Draw();
    gPad->SetLogy( true );
  }

  cv->SaveAs( Form( "Figures/RawDataTestBCO_%08i.pdf", runnumber ) );
  
  // print
  if( false )
  {
    for( const auto& [fee, bco_map]:fee_bco_map )
    {
      for( const auto& [bco, entries]:bco_map )
      { 
        std::cout 
          << "RawDataTestBCO -"
          << " fee: " << fee \
          << " bco: " << bco 
          << " entries: " << entries << std::endl; 
      }
    }
    
    std::cout << std::endl;    
  }
  
}
