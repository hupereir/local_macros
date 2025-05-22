#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <fun4all/Fun4AllUtils.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________
void CheckWaveformHistory()
{

  set_style(false);

  const std::vector<TString> fileSelection =
  {
    "/sphenix/data/data02/sphnxpro/streamhist/new_2024p002/*/HIST*.root",
    "/sphenix/data/data02/sphnxpro/streamhist/tpotonly/new_2024p002/*/HIST*.root"
  };

  // create filemanager, load file selections
  FileManager f;
  for( const auto& selection:fileSelection )
  { f.AddFiles( selection ); }

  // get single files
  const auto files = f.GetFiles();
  std::cout << "CheckWaveformHistory - files: " << files.size() << std::endl;

  struct data_t
  {
    double entries = 0;
    double mean = 0;
    double rms = 0;
  };

  std::map<int,data_t> data_map;

  int i = 0;
  for( const auto& file:files )
  {
    // std::cout << "CheckWaveformHistory - file: " << file << std::endl;

    // readout initialization
    const auto [runnumber,segment] = Fun4AllUtils::GetRunSegment(file.Data());

    // open TFile
    TFile f(file);
    auto h = static_cast<TH1*>( f.Get("h_MicromegasBCOQA_nwaveform_bco") );
    if( !h ) continue;

    data_t data;

    data.entries = h->GetEntries();
    data.mean = h->GetMean();
    data.rms = h->GetRMS();
    if( data.entries < 1 ) continue;
    // if( data.entries < 1e4 ) continue;

    std::cout << "CheckWaveformHistory - runnumber: " << runnumber << " mean: " << data.mean << " rms: " << data.rms << std::endl;
    data_map.emplace(runnumber, data);

  }


  std::cout << "CheckWaveformHistory - files: " << files.size() << " run count: " << data_map.size() << std::endl;

  auto tg_mean = new TGraphErrors();
  auto tg_rms = new TGraphErrors();
  const int color[] = { 1, 1};
  const TString label[] =
  {
    "#LT waveforms/trigger #GT",
    "RMS waveforms/trigger"
  };

  {
    int i =0;
    for( const auto& tg:{tg_mean, tg_rms} )
    {
//       tg->SetMaximum(1.05);
      tg->SetMinimum(0.);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->GetYaxis()->SetTitle( label[i] );
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
      tg_mean->SetPoint(i, runnumber, data.mean);
      tg_rms->SetPoint(i, runnumber, data.rms);
      ++i;
    }
  }

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  cv->Divide(1, 2);
  cv->cd(1);
  tg_mean->Draw("APL");

  cv->cd(2);
  tg_rms->Draw("APL");

  cv->SaveAs("Figures/WaveformHistory.pdf");
}
