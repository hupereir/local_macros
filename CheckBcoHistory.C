#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <fun4all/Fun4AllUtils.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________
void CheckBcoHistory()
{

  set_style(false);

  const std::vector<TString> fileSelection =
  {
    "/sphenix/data/data02/sphnxpro/streamhist/new_2024p002/*/HIST*.root",
    "/sphenix/data/data02/sphnxpro/streamhist/tpotonly/new_2024p002/*/HIST*.root",
    "/phenix/u/hpereira/sphenix/work/g4simulations/Rootfiles/QA/CombinedDataQA-*-0000.root"
  };

  // create filemanager, load file selections
  FileManager f;
  for( const auto& selection:fileSelection )
  { f.AddFiles( selection ); }

  // get single files
  const auto files = f.GetFiles();
  std::cout << "CheckBcoHistory - files: " << files.size() << std::endl;
//
  struct data_t
  {
    double ref = 0;
    double found_5001 = 0;
    double found_5002 = 0;
    double found = 0;
  };

  std::map<int,data_t> data_map;

  int i = 0;
  for( const auto& file:files )
  {
    std::cout << "CheckBcoHistory - file: " << file << std::endl;

    // readout initialization
    const auto [runnumber,segment] = Fun4AllUtils::GetRunSegment(file.Data());

//     if( !runnumber )
//     { std::cout << "CheckBcoHistory 0

    // open TFile
    TFile f(file);
    auto h = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_packet_stat") );
    if( !h ) continue;

    data_t data;

    data.ref = h->GetBinContent(1);
    data.found_5001 = h->GetBinContent(2);
    data.found_5002 = h->GetBinContent(3);
    data.found = h->GetBinContent(4);
    if( data.ref < 1e4 ) continue;

    std::cout << "CheckBcoHistory - runnumber: " << runnumber << " ref: " << data.ref << " found: " << data.found << " eff: " << data.found/data.ref << std::endl;

    data_map.emplace(runnumber, data);

  }


  auto tg_5001 = new TGraphErrors();
  auto tg_5002 = new TGraphErrors();
  auto tg_both = new TGraphErrors();
  const int color[] = { 1, 4, 2 };
  {
    int i =0;
    for( const auto& tg:{tg_5001, tg_5002, tg_both} )
    {
      tg->SetMaximum(1.05);
      tg->SetMinimum(0.8);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->GetYaxis()->SetTitle( "BCO fraction (data / GL1)" );
      tg->SetLineColor(color[i]);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(color[i]);

      ++i;
    }
  }

  {
    int i = 0;
    for( const auto& [runnumber, data]:data_map )
    {
      tg_5001->SetPoint(i, runnumber, data.found_5001/data.ref);
      tg_5002->SetPoint(i, runnumber, data.found_5002/data.ref);
      tg_both->SetPoint(i, runnumber, data.found/data.ref);
      ++i;
    }
  }

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

  cv->SaveAs("Figures/BcoHistory.pdf");
}
