#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <fun4all/Fun4AllUtils.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include <micromegas/MicromegasMapping.h>

#include "run_list_physics.h"

namespace
{

  // detector names
  std::vector<std::string> detector_names;

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

}

//____________________________________________
void CheckOnlineQAHistory()
{

  set_style(false);
  save_detector_names();

  // store relevant data
  class data_t
  {
    public:

    // hit charge
    double hit_charge_mean = 0;
    double hit_charge_rms = 0;

    // timed hit multiplicity
    double hit_mult_mean = 0;
    double hit_mult_rms = 0;

    // waveform (untimed) multiplicyt
    double wf_mult_mean = 0;
  };

  using data_array_t = std::array<data_t, 16>;
  using data_map_t = std::map<int, data_array_t>;
  data_map_t data_map;

  for( const auto& runnumber:RUNS::runnumbers )
  {
    const auto filename = Form( "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/TPOTMON/Run_%i-TPOTMON_0.root", runnumber );
    TFile f(filename, "READ");
    if( !f.IsOpen() ) { continue; }

    // get number of triggers
    const auto h = static_cast<TH1*>( f.Get("m_counters") );
    const double trigger_count = h->GetBinContent(3);
    if( trigger_count < 1e5 ) continue;

    auto& data_array = data_map[runnumber];

    // loop over detectors
    for( int i=0; i<16; ++i )
    {
      auto& data = data_array[i];
      const auto& detname = detector_names[i];

      // hit charge
      {
        const auto hname = Form( "m_hit_charge_%s", detname.c_str() );
        const auto h = static_cast<TH1*>( f.Get(hname) );
        if( h )
        {
          data.hit_charge_mean = h->GetMean();
          data.hit_charge_rms = h->GetRMS();
        }
      }

      // hit multiplicity
      {
        const auto hname = Form( "m_hit_multiplicity_%s", detname.c_str() );
        const auto h = static_cast<TH1*>( f.Get(hname) );
        if( h )
        {
          data.hit_mult_mean = h->GetMean();
          data.hit_mult_rms = h->GetRMS();
        }
      }

      // waveform multiplicity
      {
        const auto hname = Form( "m_wf_vs_channel_%s", detname.c_str() );
        const auto h = static_cast<TH1*>( f.Get(hname) );
        if( h )
        {
          data.wf_mult_mean = double(h->GetEntries())/trigger_count;
        }
      }
    }

  }

  // make some plots
  using tg_array_t = std::array<TGraphErrors*, 16>;
  tg_array_t tg_hit_charge_mean;
  tg_array_t tg_hit_mult_mean;
  tg_array_t tg_wf_mult_mean;
  for( int i = 0; i < 16; ++i )
  {
    {
      // hit charge
      tg_hit_charge_mean[i] = new TGraphErrors();
      tg_hit_charge_mean[i]->SetMinimum(150);
      tg_hit_charge_mean[i]->SetMaximum(350);

      tg_hit_charge_mean[i]->GetXaxis()->SetTitle( "run number" );
      tg_hit_charge_mean[i]->GetYaxis()->SetTitle( "#LTQ_{hit}#GT (ADC)" );
      tg_hit_charge_mean[i]->GetYaxis()->SetTitleOffset(1.2);
    }

    {
      // hit multiplicity
      tg_hit_mult_mean[i] = new TGraphErrors();
      tg_hit_mult_mean[i]->SetMinimum(0);
      tg_hit_mult_mean[i]->SetMaximum(2);

      tg_hit_mult_mean[i]->GetXaxis()->SetTitle( "run number" );
      tg_hit_mult_mean[i]->GetYaxis()->SetTitle( "#LTN_{hit}#GT" );
      tg_hit_mult_mean[i]->GetYaxis()->SetTitleOffset(1.2);
    }

    {
      // waveform multiplicity
      tg_wf_mult_mean[i] = new TGraphErrors();
      tg_wf_mult_mean[i]->SetMinimum(0);
      tg_wf_mult_mean[i]->SetMaximum(20);

      tg_wf_mult_mean[i]->GetXaxis()->SetTitle( "run number" );
      tg_wf_mult_mean[i]->GetYaxis()->SetTitle( "#LTN_{waveform}#GT" );
      tg_wf_mult_mean[i]->GetYaxis()->SetTitleOffset(1.2);
    }

    for( auto tg: { tg_hit_charge_mean[i], tg_hit_mult_mean[i], tg_wf_mult_mean[i] })
    {
      tg->SetLineColor(1);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(1);
      tg->SetMarkerSize(1);
    }


  }

  // loop over data
  int ipt = 0;
  for( const auto& [runnumber, data_array]:data_map )
  {
    for( int i = 0; i < 16; ++i )
    {
      tg_hit_charge_mean[i]->SetPoint(ipt, runnumber, data_array[i].hit_charge_mean);
      tg_hit_mult_mean[i]->SetPoint(ipt, runnumber, data_array[i].hit_mult_mean);
      tg_wf_mult_mean[i]->SetPoint(ipt, runnumber, data_array[i].wf_mult_mean);
    }

    ++ipt;
  }

  // make plot
  PdfDocument pdfDocument( "Figures/OnlineQAHistory.pdf" );

  {
    // hit charge
    auto cv = new TCanvas( "cv", "", 1200, 800 );
    cv->Divide(4,4);
    for( int i = 0; i < 16; ++i )
    {
      cv->cd(i+1);
      tg_hit_charge_mean[i]->Draw("APL");

      gPad->Update();
      Draw::PutText( 0.2, 0.7, detector_names[i].c_str() );
    }
    pdfDocument.Add(cv);
  }

  {
    // hit mult
    auto cv = new TCanvas( "cv1", "", 1200, 800 );
    cv->Divide(4,4);
    for( int i = 0; i < 16; ++i )
    {
      cv->cd(i+1);
      tg_hit_mult_mean[i]->Draw("APL");

      gPad->Update();
      Draw::PutText( 0.2, 0.7, detector_names[i].c_str() );
    }
    pdfDocument.Add(cv);
  }

  {
    // wf mult
    auto cv = new TCanvas( "cv2", "", 1200, 800 );
    cv->Divide(4,4);
    for( int i = 0; i < 16; ++i )
    {
      cv->cd(i+1);
      // tg_wf_mult_mean[i]->GetXaxis()->SetRangeUser(52000, 53000);
      tg_wf_mult_mean[i]->Draw("APL");

      // gPad->SetLogy(true);
      gPad->Update();
      Draw::PutText( 0.2, 0.7, detector_names[i].c_str() );
    }
    pdfDocument.Add(cv);
  }

}
