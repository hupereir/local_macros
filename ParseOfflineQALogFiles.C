#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>

#include <TPRegexp.h>

#include <fstream>

R__LOAD_LIBRARY(libRootUtilBase.so)

//_____________________________________________________-
void ParseOfflineQALogFiles()
{

  // const TString tag = "-old";
  const TString tag = "-new";
  // const TString tag;
  // const auto inputfilenames = Form( "LOG/output_OfflineQA_hp_*_full%s.txt", tag.Data() );
  // const auto inputfilenames = Form( "LOG/output_ReadCombinedData_GL1_*_full%s.txt", tag.Data() );

  // const auto inputfilenames = Form( "LOG/output_ReadCombinedData_GL1_1068_full%s.txt", tag.Data() );

  // const auto inputfilenames = Form( "LOG/output_ReadCombinedData_GL1_693_full%s.txt", tag.Data() );
  // const auto inputfilenames = Form( "LOG/output_ReadCombinedData_GL1_1154_full%s.txt", tag.Data() );
  const auto inputfilenames = Form( "LOG/output_ReadCombinedData_GL1_1403_full%s.txt", tag.Data() );

  const auto pdffilename = Form( "Figures/OfflineQA_timing%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdffilename );

  // create file manager
  FileManager fileManager( inputfilenames );

  // run number regexp
  // inputfile: LUSTRE_PHYSICS/physics//TPOT_ebdc39_physics-00049356-0000.evt
  TPRegexp runnumber_regexp("physics-(\\d+)-\\d+\\.evt");

  // timing regexp
  // EventCounter_hp::process_event - event = 109800 time (ms): 1348.86
  TPRegexp regexp( "EventCounter_hp::process_event - event = (\\d+) time \\(ms\\): (\\d+\\.\\d+)" );

  // decoder timing regex
  TPRegexp decoder_timing_regex( "SingleMicromegasPoolInput(_v2)?: per event time \\(ms\\):\\s+(\\d+\\.\\d+)" );

  // waveforms
  // SingleMicromegasPoolInput::~SingleMicromegasPoolInput - packet: 5001 wf_total: 24796148 wf_dropped_bco: 89782 wf_dropped_pool: 5898 ratio_bco: 0.0036208 ratio_pool: 0.00023786

  TPRegexp wf_5001_regex( "packet: 5001 wf_total: (\\d+)" );
  TPRegexp wf_5002_regex( "packet: 5002 wf_total: (\\d+)" );

  // decoder timing
  auto h_decoder_timing = new TH1F("h_decoder_timing", "decoder time/evt (ms)", 100, 0, 2 );
  h_decoder_timing->GetXaxis()->SetTitle( "time/RCDaq frame (ms)" );

  // mean per event time
  auto h_mean_timing = new TH1F( "h_mean_timing", "mean time per 100 events (ms) per runnumber", 100, 0, 500 );
  h_mean_timing->GetXaxis()->SetTitle( "mean time/100 events per run number (ms)" );

  // mean number of waveform per packet
  auto h_waveforms_5001 = new TH1F( "h_wf_5001", "number of waveforms/event (packet 5001)", 100, 0, 200 );
  h_waveforms_5001->GetXaxis()->SetTitle( "number of waveforms/event" );

  auto h_waveforms_5002 = new TH1F( "h_wf_5002", "number of waveforms/event (packet 5002)", 100, 0, 200 );
  h_waveforms_5002->GetXaxis()->SetTitle( "number of waveforms/event" );

  // loop over selected files
  const auto filenames = fileManager.GetFiles();
  for( const auto& filename:filenames )
  {
    std::cout << "ParseOfflineQALogFiles - processing " << filename << std::endl;

    int runnumber = 0;
    TH2* h_time = nullptr;

    std::ifstream in( filename.Data() );
    std::string line;
    while( std::getline(in,line) )
    {
      if( !h_time )
      {
        // run number
        auto matches = runnumber_regexp.MatchS( line.c_str() );
        if( matches->GetLast() >= 1 )
        {
          runnumber =  static_cast<TObjString*>(matches->At(1))->GetString().Atof();
          std::cout << "runnumber: " << runnumber << std::endl;
          const auto hname = Form( "htime_%i", runnumber );
          const auto title = Form( "time (ms)/100 events - run number: %i", runnumber );
          h_time = new TH2F( hname, title, 100, 0, 500000, 100, 0, 500 );
          h_time->SetStats(0);
          h_time->GetXaxis()->SetTitle( "event" );
          h_time->GetYaxis()->SetTitle( "time/100 ev (ms)" );
          h_time->GetYaxis()->SetTitleOffset( 1.5 );
          continue;
        }
      }

      // timing
      if( h_time )
      {
        auto matches = regexp.MatchS( line.c_str() );
        if( matches->GetLast() >= 2 )
        {
          const int event =  static_cast<TObjString*>(matches->At(1))->GetString().Atof();
          const float time =  static_cast<TObjString*>(matches->At(2))->GetString().Atof();
          h_time->Fill(event,time);

          continue;
        }
      }

      // decoder timing
      {
        auto matches = decoder_timing_regex.MatchS( line.c_str() );
        if( matches->GetLast() >= 2 )
        {
          const float time =  static_cast<TObjString*>(matches->At(2))->GetString().Atof();
          h_decoder_timing->Fill(time);
          std::cout << "Decoder timing: " << time << std::endl;
        }
      }

      // number of waveforms
      {
        auto matches = wf_5001_regex.MatchS( line.c_str() );
        if( matches->GetLast() >= 1 )
        {
          const float wf =  static_cast<TObjString*>(matches->At(1))->GetString().Atof()/5e5;
          h_waveforms_5001->Fill(wf);
          std::cout << "wf: " << wf << std::endl;
        }
      }

      {
        auto matches = wf_5002_regex.MatchS( line.c_str() );
        if( matches->GetLast() >= 1 )
        {
          const float wf =  static_cast<TObjString*>(matches->At(1))->GetString().Atof()/5e5;
          h_waveforms_5002->Fill(wf);
        }
      }

    }
    if( h_time && h_time->GetEntries()>100 )
    {
      h_mean_timing->Fill( h_time->GetMean(2) );
    }

    if( true && h_time && h_time->GetEntries()>100 )
    {
      auto cv_name = Form( "cv_%i", runnumber );
      auto cv = new TCanvas( cv_name, cv_name, 1200, 900 );
      h_time->Draw();
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.1);
      pdfDocument.Add(cv);
    }
  }

  if( true && h_mean_timing->GetEntries()>0  )
  {
    auto cv_name = TString( "cv_mean_timing" );
    auto cv = new TCanvas( cv_name, cv_name, 1200, 900 );
    h_mean_timing->Draw();
    pdfDocument.Add(cv);
  }


  if( true && h_decoder_timing->GetEntries()>0  )
  {
    auto cv_name = TString( "cv_decoder_timing" );
    auto cv = new TCanvas( cv_name, cv_name, 1200, 900 );
    h_decoder_timing->Draw();
    pdfDocument.Add(cv);
  }

  if( true && h_waveforms_5001->GetEntries()>0  )
  {
    auto cv_name = TString( "cv_wf" );
    auto cv = new TCanvas( cv_name, cv_name, 1200, 900 );
    cv->Divide(2,1);
    cv->cd(1);
    h_waveforms_5001->Draw();

    cv->cd(2);
    h_waveforms_5002->Draw();
    pdfDocument.Add(cv);
  }

}
