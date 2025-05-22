#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <fun4all/Fun4AllUtils.h>

// #include "run_list_physics.h"
#include "run_list_test.h"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________
void CheckBcoHistory_QA()
{

  set_style(false);

  struct data_t
  {
    // BCO matching statistics
    double ref = 0;
    double found_5001 = 0;
    double found_5002 = 0;
    double found = 0;

    // dropped waveforms statistics
    double wf_ref_5001 = 0;
    double wf_ref_5002 = 0;
    double wf_dropped_bco_5001 = 0;
    double wf_dropped_bco_5002 = 0;

    double wf_dropped_pool_5001 = 0;
    double wf_dropped_pool_5002 = 0;

  };

  std::map<int,data_t> data_map;

  int total_runs = 0;
  int good_runs = 0;

  int i = 0;

  const TString tag = "-new";

  // create plots
  PdfDocument pdfDocument(Form("Figures/BcoHistory_QA%s.pdf", tag.Data()));

  for( const auto& runnumber:RUNS::runnumbers )
  {

    // path
    const int range_min = (int)(runnumber/100)*100;
    const int range_max = range_min+100;
//     const auto path = Form( "/sphenix/data/data02/sphnxpro/testbed/single_streamhist/ana450_2024p008/run_%08i_%08i/", range_min, range_max );
//     const auto filename = Form( "%s/HIST_DST_STREAMING_EVENT_TPOT_run2pp_ana450_2024p008-%08i-00000.root", path, runnumber );

    // const auto filename = Form( "Rootfiles/QA/MicromegasClusterQA_GL1%s/MicromegasClusterQA-%08i-0000.root", tag.Data(), runnumber );
    const auto filename = Form( "Rootfiles/QA/CombinedData_GL1%s/CombinedDataQA-%08i-0000.root", tag.Data(), runnumber );

    std::cout << "CheckBcoHistory_QA -"
      << " runnumber: " << runnumber
      << " filename: " << filename
      << std::endl;

    // open TFile
    TFile f(filename, "READ");
    if( !f.IsOpen() ) { continue; }

    data_t data;

    {
      // BCO statistics
      auto h = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_packet_stat") );
      if( !h ) continue;

      data.ref = h->GetBinContent(1);
      data.found_5001 = h->GetBinContent(2);
      data.found_5002 = h->GetBinContent(3);
      data.found = h->GetBinContent(4);

      // skip if too few events
      // if( data.ref < 5e5 ) continue;
      if( data.ref < 1e5 ) continue;
      std::cout << "CheckBcoHistory - runnumber: " << runnumber << " ref: " << data.ref << " found: " << data.found << " eff: " << data.found/data.ref << std::endl;
    }

    {
      // dropped waveforms statistics
      auto h_ref = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_waveform_count_total") );
      data.wf_ref_5001 = h_ref->GetBinContent(1);
      data.wf_ref_5002 = h_ref->GetBinContent(2);

      // auto h_dropped_bco = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_waveform_count_dropped") );
      auto h_dropped_bco = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_waveform_count_dropped_bco") );
      if( h_dropped_bco )
      {
        data.wf_dropped_bco_5001 = h_dropped_bco->GetBinContent(1);
        data.wf_dropped_bco_5002 = h_dropped_bco->GetBinContent(2);
      }

      auto h_dropped_pool = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_waveform_count_dropped_pool") );
      if( h_dropped_pool )
      {
        data.wf_dropped_pool_5001 = h_dropped_pool->GetBinContent(1);
        data.wf_dropped_pool_5002 = h_dropped_pool->GetBinContent(2);
      }
    }

    // insert in map
    data_map.emplace(runnumber, data);

  }


  // dropped BCO fraction
  auto tg_5001 = new TGraphErrors();
  auto tg_5002 = new TGraphErrors();
  auto tg_both = new TGraphErrors();
  const int color[] = { 1, 4, 2 };
  {
    int i =0;
    for( const auto& tg:{tg_5001, tg_5002, tg_both} )
    {
      tg->SetMaximum(1.01);
      tg->SetMinimum(0.9);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->GetYaxis()->SetTitle( "BCO fraction (data / GL1)" );
      tg->SetLineColor(color[i]);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(color[i]);

      ++i;
    }
  }

  // dropped waveform fraction (bco)
  auto tg_wf_dropped_bco_5001 = new TGraphErrors();
  auto tg_wf_dropped_bco_5002 = new TGraphErrors();
  {
    int i =0;
    for( const auto& tg:{tg_wf_dropped_bco_5001, tg_wf_dropped_bco_5002} )
    {
      tg->SetMinimum(0);
      tg->SetMaximum(1.);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->GetYaxis()->SetTitle( "dropped waveform fraction (bco)" );
      tg->SetLineColor(color[i]);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(color[i]);

      ++i;
    }
  }

  // dropped waveform fraction (bco)
  auto tg_wf_dropped_pool_5001 = new TGraphErrors();
  auto tg_wf_dropped_pool_5002 = new TGraphErrors();
  {
    int i =0;
    for( const auto& tg:{tg_wf_dropped_pool_5001, tg_wf_dropped_pool_5002} )
    {
      tg->SetMinimum(0);
      tg->SetMaximum(1.);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->GetYaxis()->SetTitle( "dropped waveform fraction (pool)" );
      tg->SetLineColor(color[i]);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(color[i]);

      ++i;
    }
  }

  /// fill all tgraphs
  {
    int i = 0;
    std::array<int,2> i_wf = {};
    for( const auto& [runnumber, data]:data_map )
    {

      tg_5001->SetPoint(i, runnumber, data.found_5001/data.ref);
      tg_5002->SetPoint(i, runnumber, data.found_5002/data.ref);
      tg_both->SetPoint(i, runnumber, data.found/data.ref);
      ++i;
      ++total_runs;
      if( data.found/data.ref>0.98 )
      {
        ++good_runs;
      } else {
        std::cout << "CheckBcoHistory_QA - bad run: " << runnumber << std::endl;
      }

      if( data.wf_ref_5001>0 )
      {
        tg_wf_dropped_bco_5001->SetPoint(i_wf[0], runnumber, data.wf_dropped_bco_5001/data.wf_ref_5001 );
        tg_wf_dropped_pool_5001->SetPoint(i_wf[0], runnumber, data.wf_dropped_pool_5001/data.wf_ref_5001 );
        ++i_wf[0];
      }

      if( data.wf_ref_5002>0 )
      {
        tg_wf_dropped_bco_5002->SetPoint(i_wf[1], runnumber, data.wf_dropped_bco_5002/data.wf_ref_5002 );
        tg_wf_dropped_pool_5002->SetPoint(i_wf[1], runnumber, data.wf_dropped_pool_5002/data.wf_ref_5002 );
        ++i_wf[1];
      }

      std::cout << "CheckBcoHistory_QA -"
        << " runnumber: " << runnumber
        << " ref_5001: " << data.wf_ref_5001
        << " dropped_5001: " << data.wf_dropped_bco_5001/data.wf_ref_5001
        << " dropped_5001_pool: " << data.wf_dropped_pool_5001/data.wf_ref_5001
        << std::endl;

      std::cout << "CheckBcoHistory_QA -"
        << " runnumber: " << runnumber
        << " ref_5002: " << data.wf_ref_5002
        << " dropped_5002: " << data.wf_dropped_bco_5002/data.wf_ref_5002
        << " dropped_pool_5002: " << data.wf_dropped_pool_5002/data.wf_ref_5002
        << std::endl;
    }
  }

  std::cout << "total_runs: " << total_runs << " good_runs: " << good_runs << std::endl;

  {
    // dropped packets
    auto cv = new TCanvas( "cv", "cv", 1200, 800 );
    tg_5001->Draw("APL");
    tg_5002->Draw("PL");
    tg_both->Draw("PL");

    auto legend = new TLegend( 0.15, 0.24, 0.82, 0.39, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    legend->AddEntry(tg_5001, "pkt 5001", "L" );
    legend->AddEntry(tg_5002, "pkt 5002", "L" );
    legend->AddEntry(tg_both, "all", "L" );
    legend->Draw();
    pdfDocument.Add(cv);
  }

  {
    // dropped waveforms (bco)
    auto cv = new TCanvas( "cv2", "cv2", 1200, 800 );
    tg_wf_dropped_bco_5001->Draw("APL");
    tg_wf_dropped_bco_5002->Draw("PL");

    auto legend = new TLegend( 0.15, 0.76, 0.42, 0.9, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(0);

    legend->AddEntry(tg_wf_dropped_bco_5001, "pkt 5001", "L" );
    legend->AddEntry(tg_wf_dropped_bco_5002, "pkt 5002", "L" );
    legend->Draw();

    if( true )
    {
      auto line = Draw::VerticalLine(gPad, 51209);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }

    if( false )
    {
      auto line = Draw::VerticalLine(gPad, 51865);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }

    if( true )
    {
      auto line = Draw::VerticalLine(gPad, 52077);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }
    pdfDocument.Add(cv);
  }


  {
    // dropped waveforms (pool)
    auto cv = new TCanvas( "cv3", "cv3", 1200, 800 );
    tg_wf_dropped_pool_5001->Draw("APL");
    tg_wf_dropped_pool_5002->Draw("PL");

    auto legend = new TLegend( 0.15, 0.76, 0.42, 0.9, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(0);

    legend->AddEntry(tg_wf_dropped_pool_5001, "pkt 5001", "L" );
    legend->AddEntry(tg_wf_dropped_pool_5002, "pkt 5002", "L" );
    legend->Draw();

    if( true )
    {
      auto line = Draw::VerticalLine(gPad, 51209);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }

    if( false )
    {
      auto line = Draw::VerticalLine(gPad, 51865);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }

    if( true )
    {
      auto line = Draw::VerticalLine(gPad, 52077);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw();
    }
    pdfDocument.Add(cv);
  }

}
