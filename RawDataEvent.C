#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

namespace
{

  // test if a given set contains a given value
  template<class T>
    bool contains( const std::set<T>& container, const T& value )
  { return container.find( value ) != container.end(); }
  
  // detector names
  std::vector<std::string> detector_names;

  // event processing
  void process_event( PdfDocument& pdfDocument, uint64_t lvl1_bco, const MicromegasRawDataEvaluation::Waveform::List& waveforms )
  {
//     using hist_pointer = std::unique_ptr<TH1>;
//     std::array<hist_pointer,16> histograms;
    std::array<TH1*,16> histograms;
        
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& detector_name = detector_names[detid];
      
      // create histogram
      const auto hname = Form( "h_channel_%lu_%s", lvl1_bco, detector_name.c_str() );
      const auto htitle = Form( "BCO: %lu %s",  lvl1_bco, detector_name.c_str() );
      auto h = new TH1D( hname, htitle, 256, 0, 256 );
      h->GetXaxis()->SetTitle( "strip" );
      h->GetYaxis()->SetTitle( "adc" );
      h->SetMaximum(1200);

      // histograms[detid].reset(h);
      histograms[detid] = h;
      
    }
                
    // loop over waveforms and fill
    for( const auto& wf:waveforms )
    { 
      if( wf.layer == 0 ) continue;
      
      const auto& itile = wf.tile;
      const auto ilayer = wf.layer-55;
      const int detid = itile + 8*ilayer;
      histograms[detid]->Fill( wf.strip, wf.adc_max );
    }

    // create canvas
    auto cvname = Form( "cv_%lu",  lvl1_bco );
    // auto cv = std::make_unique<TCanvas>( cvname, cvname, 900, 500 );
    auto cv = new TCanvas( cvname, cvname, 1200, 900 );
    cv->Divide( 4,4 );
    
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& h = histograms[detid];
      cv->cd(detid+1);
      h->Draw( "HIST" );    
    }

    // pdfDocument.Add(cv.get());    
    pdfDocument.Add(cv);    
  }
}

//_____________________________________________________________________________
TString RawDataEvent( int runNumber = 20445 )
{
  

  const std::set<uint64_t> selected_events;
  
  // get detector names that match tile/layer
  MicromegasMapping mapping;
  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
  {
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
    const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
    detector_names.push_back( std::move( name ) );
  }

  const TString inputFile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataEvent-%08i-0000.pdf", runNumber );
  
  std::cout << "RawDataEvaluation - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataEvaluation - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // keep track of current waveforms
  MicromegasRawDataEvaluation::Waveform::List waveforms;
  
  // keep track of all waveforms associated to a given bco 
  using waveform_map_t = std::map<uint64_t, MicromegasRawDataEvaluation::Waveform::List>;
  waveform_map_t waveform_map;
  
  // loop over tree entries
  for( int i = 0; i < 500; ++i )
  // for( int i = 0; i < tree->GetEntries(); ++i )
  {
    if( !(i%100) ) std::cout << "RawDataEvent - entry: " << i << std::endl;
    tree->GetEntry(i);
    
    // loop over waveforms
    for( const auto& waveform:container->waveforms )
    { waveform_map[waveform.lvl1_bco].push_back( waveform ); }    
  }

  for( const auto& [lvl1_bco,waveforms]:waveform_map)
  {
    std::cout
      << "RawDataEvent -"
      << " bco: " << lvl1_bco 
      << " waveforms: " << waveforms.size() 
      << std::endl;
        
    // make plot
    process_event( pdfDocument, lvl1_bco, waveforms );          
  }
  
  return pdfFile;
}
