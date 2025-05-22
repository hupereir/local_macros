#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <fun4all/Fun4AllUtils.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{

  // streamer for lists
  template <class T>
  std::ostream& operator<<(std::ostream& o, const std::set<T>& list)
  {
    if (list.empty())
    {
      o << "{}";
    }
    else
    {

      int count = 0;

      const bool is_hex = (o.flags() & std::ios_base::hex);
      o << "{ " << std::endl << " ";
      bool first = true;
      for (const auto& value : list)
      {
        if (count>0)
        {
          o << ", ";
        }

        if (is_hex)
        {
          o << "0x";
        }
        o << value;

        if( ++count == 10 )
        {
          count = 0;
          o << ", " << std::endl << " ";
        }

        first = false;
      }
      o << std::endl << " }";
    }
    return o;
  }


}

//____________________________________________
void CheckRunHistory()
{

  set_style(false);

  // get all processed runs
  std::set<int> runnumbers_data;
  {
    const std::vector<TString> fileSelection =
    {
      "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/run_*/*-00000.root",
      "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/fast/run_*/*-00000.root",
      "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002/run_*/*-00000.root"
    };

    // create filemanager, load file selections
    FileManager f;
    for( const auto& selection:fileSelection )
    { f.AddFiles( selection ); }

    // get single files
    const auto files = f.GetFiles();
    for( const auto& file:files )
    {
      const auto [runnumber,segment] = Fun4AllUtils::GetRunSegment(file.Data());
      if(!runnumber)
      {
        std::cout << "CheckRunHistory - could not get run number: " << file << std::endl;
        continue;
      }

      // todo: check file size or event numbers to cut out small runs.
      {
        TFile f(file);
        if( !f.IsOpen() ) continue;
        auto tree = static_cast<TTree*>(f.Get("T"));
        if( !tree ) continue;

        if( tree->GetEntries()<100) { continue; }

        std::cout << "CheckRunHistory - runnumber: " << runnumber << " entries: " << tree->GetEntries() << std::endl;
      }
      runnumbers_data.insert(runnumber);
    }

    std::cout << "CheckRunHistory - runnumbers_data: " << runnumbers_data.size() << std::endl;

  }

  std::cout << "std::vector<int> runnumbers = " << runnumbers_data << "; " << std::endl;

  // get all runs for which there are histograms.
  std::set<int> runnumbers_hist;
  {
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
    for( const auto& file:files )
    {
      const auto [runnumber,segment] = Fun4AllUtils::GetRunSegment(file.Data());
      if( runnumber ) { runnumbers_hist.insert(runnumber); }
      else { std::cout << "CheckRunHistory - could not get run number: " << file << std::endl; }
    }

    std::cout << "CheckRunHistory - runnumbers_hist: " << runnumbers_hist.size() << std::endl;

  }

  // store
  auto tg_data = new TGraphErrors();
  auto tg_hist = new TGraphErrors();
  const int color[] = { 1, 2 };
  {
    int i =0;
    for( const auto& tg:{tg_data, tg_hist} )
    {
      tg->SetMaximum(1.05);
      tg->SetMinimum(0);
      tg->GetXaxis()->SetTitle( "run number" );
      tg->SetLineColor(color[i]);
      tg->SetLineWidth(2);

      tg->SetMarkerColor(color[i]);
      tg->SetMarkerStyle(20);

      ++i;
    }
  }

  {
    int i = 0;
    for( const auto& runnumber:runnumbers_data )
    { tg_data->SetPoint( i++, runnumber ,1 ); }
  }

  {
    int i = 0;
    for( const auto& runnumber:runnumbers_hist )
    { tg_hist->SetPoint( i++, runnumber ,1 ); }
  }

  auto cv = new TCanvas( "cv", "cv", 1200, 800 );
  tg_data->Draw("AP");
  tg_hist->Draw("P");

  auto legend = new TLegend( 0.15, 0.24, 0.82, 0.39, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(tg_data, "data  available", "P" );
  legend->AddEntry(tg_hist, "QA histogram  available", "P" );
  legend->Draw();

  cv->SaveAs("Figures/RunHistory.pdf");
}
