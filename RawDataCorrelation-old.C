#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>
#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

namespace
{
  
  // detector names
  std::vector<std::string> detector_names;
  
  std::array<TH2*,16> histograms;

  std::array<TH2*,2> histograms_layer;

  MicromegasMapping mapping;

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

    for( int ilayer = 0; ilayer <2; ++ilayer )
    {
      {
        const auto hname = Form( "h_channel_%i", ilayer+55 );
        const auto htitle = Form( "correlation %s", (ilayer==0) ? "#phi":"z" );
        auto h = new TH2I( hname, htitle, 256, 0, 256, 256, 0, 256 );
        h->GetXaxis()->SetTitle( "strip" );
        h->GetYaxis()->SetTitle( "strip" );
        histograms_layer[ilayer] = h;
      }
      
      for( int itile = 0; itile <8; ++itile )
      {
        const int detid = itile + 8*ilayer;
        const auto& detector_name = detector_names[detid];
        
        // create histogram
        const auto hname = Form( "h_channel_%s", detector_name.c_str() );
        const auto htitle = Form( "correlation %s", detector_name.c_str() );
        auto h = new TH2I( hname, htitle, 256, 0, 256, 256, 0, 256 );
        h->GetXaxis()->SetTitle( "strip" );
        h->GetYaxis()->SetTitle( "strip" );
        histograms[detid] = h;
      }
    }
  }
  
  // event processing
  void process_event( uint64_t /*lvl1_bco*/, const MicromegasRawDataEvaluation::Waveform::List& waveforms )
  {

    for( size_t i = 0; i < waveforms.size(); ++i )
    {
      const auto& first = waveforms[i];
      if( !first.is_signal ) continue;
      if( first.layer == 0 ) continue;

      const int detid_first = first.tile + 8*(first.layer-55 );

      const auto first_strip =first.strip;
      // const auto first_strip = mapping.get_physical_strip( first.fee_id, first.channel );
      
      for( size_t j=i+1; j < waveforms.size(); ++j )
      {
        const auto& second = waveforms[j];
        if( !second.is_signal ) continue;
        if( second.layer == 0 ) continue;
        
        const int detid_second = second.tile + 8*(second.layer-55 );
        if( detid_first != detid_second ) continue;

        const auto second_strip = second.strip;
        // const auto second_strip = mapping.get_physical_strip( second.fee_id, second.channel );
        
        {
          auto h = histograms[detid_first];
          h->Fill( first_strip, second_strip );
          h->Fill( second_strip, first_strip );
        }
        
        {
          auto h = histograms_layer[first.layer-55];
          h->Fill( first_strip, second_strip );
          h->Fill( second_strip, first_strip );
        }

      }
    }   
  }
  
  void save_histograms( PdfDocument& pdfDocument )
  {
    
    for( const auto& h:histograms_layer )
    {
      auto cvname = Form( "cv_%s", h->GetName());
      auto cv = new TCanvas( cvname, cvname, 900,500 );
      cv->SetRightMargin(0.2);
      h->Draw( "colz" );      
      pdfDocument.Add( cv );
    }

    for( const auto& h:histograms )
    {
      auto cvname = Form( "cv_%s", h->GetName());
      auto cv = new TCanvas( cvname, cvname, 900,500 );
      cv->SetRightMargin(0.2);
      h->Draw( "colz" );
      gPad->SetLogz(true);
      pdfDocument.Add( cv );
    }
  }
  
}

//_____________________________________________________________________________
TString RawDataCorrelation()
{
  
  // TODO: should get the bcoid directly from the rootfile
  
//   int runNumber = 9443;
//   int runNumber = 9467;
  int runNumber = 14094;  

  const TString inputFile = Form( "DST/MicromegasRawDataEvaluation-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataCorrelation-%08i-0000.pdf", runNumber );
  
  std::cout << "RawDataEvaluation - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataEvaluation - pdfFile: " << pdfFile << std::endl;

  save_detector_names();
  create_histograms();
  
  PdfDocument pdfDocument( pdfFile );
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // keep track of current waveforms
  MicromegasRawDataEvaluation::Waveform::List waveforms;

  // keep track of current bco
  uint64_t lvl1_bco = 0;
  
  // loop over tree entries
  for( int i = 0; i < tree->GetEntries(); ++i )
  // for( int i = 0; i < 500; ++i )
  {
    tree->GetEntry(i);
    
    // loop over waveforms
    for( const auto& waveform:container->waveforms )
    {
      // check bco
      if( std::abs<int64_t>(waveform.lvl1_bco-lvl1_bco)>1 )
      {
        std::cout << "RawDataCorrelation -"
          << " bco: " << lvl1_bco 
          << " current: " << waveform.lvl1_bco 
          << " difference: " << std::abs<int64_t>(waveform.lvl1_bco-lvl1_bco)
          << " waveforms: " << waveforms.size() 
          << std::endl;
        
        // this is a new event
        if( !waveforms.empty() ) 
        {
          
          // std::cout << "bco: " << lvl1_bco << " waveforms: " << waveforms.size() << std::endl;

          // make plot
          process_event( lvl1_bco, waveforms );
          
        }
        
        // clear
        waveforms.clear();
        lvl1_bco = waveform.lvl1_bco;
      }
      
      // store waveform
      waveforms.push_back( waveform );
    }
    
  }
  
  if( !waveforms.empty() ) 
  { process_event( lvl1_bco, waveforms ); }

  // save histograms
  save_histograms( pdfDocument );
  
  std::cout << "done." << std::endl;
  
  return pdfFile;
}
