#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
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
      // std::cout << "save_detector_names - " << detector_names.size() << " " << name << std::endl;
      detector_names.push_back( std::move( name ) );
    }
  }
}

//____________________________________________
void CheckRawDataEfficiency_QA()
{

  set_style(false);
  save_detector_names();

  // store relevant data
  class data_t
  {
    public:

    double clusters_ref = 0;
    double clusters_found = 0;
  };

  using data_array_t = std::array<data_t, 16>;
  using data_map_t = std::map<int, data_array_t>;
  data_map_t data_map;

  // const TString tag = "";
  const TString tag = "-new2";
  const auto pdfname = Form("Figures/RawDataEfficiency_QA%s.pdf", tag.Data());
  const auto rootfilename = Form("Rootfiles/RawDataEfficiency_QA%s.root", tag.Data());

  std::cout << "CheckRawDataEfficiency_QA - pdfname: " << pdfname << std::endl;
  std::cout << "CheckRawDataEfficiency_QA - rootfilename: " << rootfilename << std::endl;

  PdfDocument pdfDocument(pdfname);
  RootFile rootFile( rootfilename );

  for( const auto& runnumber:RUNS::runnumbers )
  {
    const auto filename = Form("Rootfiles/QA/MicromegasClusterQA_GL1%s/MicromegasClusterQA-%08i-0000.root", tag.Data(), runnumber);
    TFile f(filename, "READ");
    if( !f.IsOpen() ) { continue; }

    // load histograms
    const auto href = static_cast<TH1*>(f.Get("h_MicromegasClusterQA_clustercount_ref"));
    const auto hfound = static_cast<TH1*>(f.Get("h_MicromegasClusterQA_clustercount_found"));
    if( !(href&&hfound)) { continue; }

    // insert data in map
    auto& data_array = data_map[runnumber];

    // loop over detectors
    for( int i=0; i<16; ++i )
    {
      auto& data = data_array[i];
      data.clusters_ref = href->GetBinContent(i+1);
      data.clusters_found = hfound->GetBinContent(i+1);
    }
  }

  // create per detector TGraphErrors
  using tg_array_t = std::array<TGraphErrors*, 16>;
  tg_array_t tg_eff;
  for( int i = 0; i < 16; ++i )
  {
    tg_eff[i] = new TGraphErrors();
    tg_eff[i]->SetName( Form("tg_eff_%i", i));
    tg_eff[i]->SetMinimum(0);
    tg_eff[i]->SetMaximum(1.1);

    tg_eff[i]->GetXaxis()->SetTitle( "run number" );
    tg_eff[i]->GetYaxis()->SetTitle( "efficiency" );
    tg_eff[i]->GetYaxis()->SetTitleOffset(1);

    tg_eff[i]->SetLineColor(1);
    tg_eff[i]->SetLineWidth(2);
    tg_eff[i]->SetMarkerStyle(20);
    tg_eff[i]->SetMarkerColor(1);
    tg_eff[i]->SetMarkerSize(1);

    rootFile.Add(tg_eff[i]);
  }

  // loop over data and fill TGraphs
  std::array<int,16> ipt={};
  for( const auto& [runnumber, data_array]:data_map )
  {
    for( int i = 0; i < 16; ++i )
    {

//       if( i==0 )
//       { std::cout << "runnumber: " << runnumber << " ref: " << data_array[i].clusters_ref << std::endl; }

      if( data_array[i].clusters_ref>500 )
      {
        const double eff = data_array[i].clusters_found/data_array[i].clusters_ref;
        const double err = std::sqrt(eff*(1.0-eff)/data_array[i].clusters_found);
        tg_eff[i]->SetPoint(ipt[i], runnumber, eff );
        tg_eff[i]->SetPointError(ipt[i], 0, err);
        ++ipt[i];
      }
    }
  }

  for( int i = 0; i < 16; ++i )
  { std::cout << "detector: " << detector_names[i] << " entries: " << ipt[i] << std::endl; }

  if( true )
  {

    // draw horizontal line for scan ZS efficiencies
    constexpr std::array<double,16> eff = {{0.617211, 0.597985, 0.553313, 0.676564, 0.593103, 0.638976, 0.534088, 0.653811, 0.406694, 0.537138, 0.621489, 0.612882, 0.513639, 0.576252, 0.572684, 0.582114}};
    constexpr std::array<double,16> err = {{0.001362, 0.001029, 0.001068, 0.001107, 0.001157, 0.001109, 0.001158, 0.001020, 0.001272, 0.001322, 0.001424, 0.001272, 0.001350, 0.001305, 0.001484, 0.001355}};

    // raw data efficiency
    auto cv = new TCanvas( "cv", "", 1200, 800 );
    cv->Divide(4,4);
    for( int i = 0; i < 16; ++i )
    {
      cv->cd(i+1);
      gPad->SetTopMargin( 0.01 );
      gPad->SetLeftMargin( 0.1 );
      gPad->SetRightMargin( 0.01 );

      tg_eff[i]->GetXaxis()->SetRange();
      tg_eff[i]->Draw("APL");

      gPad->Update();
      if( true )
      {
        auto line = Draw::VerticalLine(gPad, 51209);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        // auto line = Draw::VerticalLine(gPad, 51865);
        // auto line = Draw::VerticalLine(gPad, 51881);
        auto line = Draw::VerticalLine(gPad, 52077);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( false )
      {
        auto line = Draw::VerticalLine(gPad, 52361);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        auto line = Draw::HorizontalLine(gPad, eff[i]);
        line->SetLineColor(4);
        line->SetLineWidth(2);
        line->Draw();
      }

      Draw::PutText( 0.15, 0.85, detector_names[i].c_str(), 0.1 );
    }

    pdfDocument.Add(cv);
  }

  if( false )
  {

    // raw data efficiency
    auto cv = new TCanvas( "cv", "", 1200, 800 );
    cv->Divide(4,4);
    for( int i = 0; i < 16; ++i )
    {
      cv->cd(i+1);
      gPad->SetTopMargin( 0.01 );
      gPad->SetLeftMargin( 0.1 );
      gPad->SetRightMargin( 0.01 );
      tg_eff[i]->GetXaxis()->SetRangeUser( 52330, 52380 );
      tg_eff[i]->Draw("APL");

      gPad->Update();
      if( true )
      {
        auto line = Draw::VerticalLine(gPad, 52348);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        auto line = Draw::VerticalLine(gPad, 52361);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }


      Draw::PutText( 0.15, 0.85, detector_names[i].c_str(), 0.1 );
    }

    pdfDocument.Add(cv);
  }

  // NCOP amd NCOZ
  if( true )
  {

    // draw horizontal line for scan ZS efficiencies
    constexpr std::array<double,16> eff_zs = {{0.617211, 0.597985, 0.553313, 0.676564, 0.593103, 0.638976, 0.534088, 0.653811, 0.406694, 0.537138, 0.621489, 0.612882, 0.513639, 0.576252, 0.572684, 0.582114}};
    constexpr std::array<double,16> err_zs = {{0.001362, 0.001029, 0.001068, 0.001107, 0.001157, 0.001109, 0.001158, 0.001020, 0.001272, 0.001322, 0.001424, 0.001272, 0.001350, 0.001305, 0.001484, 0.001355}};

    constexpr std::array<double,16> eff_nzs = {{0.000000, 0.000000, 0.000000, 0.749754, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.699109, 0.000000, 0.000000, 0.000000, 0.000000}};
    constexpr std::array<double,16> err_nzs = {{0.000000, 0.000000, 0.000000, 0.001981, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.002201, 0.000000, 0.000000, 0.000000, 0.000000}};

    // raw data efficiency
    auto cv = new TCanvas( "cv", "", 800, 800 );
    cv->Divide(1,2);

    int icv = 0;
    for( int i:{3,11} )
    {
      cv->cd(++icv);

      gPad->SetTopMargin( 0.01 );
      gPad->SetLeftMargin( 0.1 );
      gPad->SetRightMargin( 0.01 );

      tg_eff[i]->GetXaxis()->SetRange();
      tg_eff[i]->Draw("APL");

      gPad->Update();
      if( true )
      {
        auto line = Draw::VerticalLine(gPad, 51209);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        // auto line = Draw::VerticalLine(gPad, 51865);
        // auto line = Draw::VerticalLine(gPad, 51881);
        auto line = Draw::VerticalLine(gPad, 52077);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        auto line = Draw::HorizontalLine(gPad, eff_zs[i]);
        line->SetLineColor(4);
        line->SetLineWidth(2);
        line->Draw();
      }

      if( true )
      {
        auto line = Draw::HorizontalLine(gPad, eff_nzs[i]);
        line->SetLineColor(2);
        line->SetLineWidth(2);
        line->Draw();
      }

      Draw::PutText( 0.15, 0.85, detector_names[i].c_str(), 0.1 );

    }

    pdfDocument.Add(cv);

  }

}
