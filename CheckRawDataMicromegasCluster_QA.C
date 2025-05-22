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
      // std::cout << "save_detector_names - " << detector_names.size() << " " << name << std::endl;
      detector_names.push_back( std::move( name ) );
    }
  }
}

//____________________________________________
void CheckRawDataMicromegasCluster_QA()
{

  set_style(false);
  save_detector_names();

  // store relevant data
  class data_t
  {
    public:

    double entries = 0;
    double cluster_multiplicity_mean = 0;
    double cluster_multiplicity_rms = 0;

    double cluster_size_mean = 0;
    double cluster_size_rms = 0;

    double cluster_charge_mean = 0;
    double cluster_charge_rms = 0;
  };

  using data_array_t = std::array<data_t, 16>;
  using data_map_t = std::map<int, data_array_t>;
  data_map_t data_map;

  PdfDocument pdfDocument( "Figures/RawDataMicromegasCluster_QA-new.pdf" );

  for( const auto& runnumber:RUNS::runnumbers )
  {
    // const auto filename = Form( "Rootfiles/QA/MicromegasClusterQA_GL1/MicromegasClusterQA-%08i-0000.root", runnumber );
    const auto filename = Form( "Rootfiles/QA/MicromegasClusterQA_GL1-new/MicromegasClusterQA-%08i-0000.root", runnumber );
    TFile f(filename, "READ");
    if( !f.IsOpen() ) { continue; }

    // std::cout << "filename: " << filename << std::endl;

    const auto h_mult = static_cast<TH2*>( f.Get("h_MicromegasClusterQA_cluster_multiplicity" ) );
    const auto h_size = static_cast<TH2*>( f.Get("h_MicromegasClusterQA_cluster_size" ) );
    const auto h_charge = static_cast<TH2*>( f.Get("h_MicromegasClusterQA_cluster_charge" ) );

    auto& data_array = data_map[runnumber];

    // loop over detectors
    for( int i=0; i<16; ++i )
    {
      auto& data = data_array[i];
      if(h_mult)
      {
        auto h = h_mult->ProjectionY( "_proj", i+1, i+1 );
        data.entries = h->Integral();
        data.cluster_multiplicity_mean = h->GetMean();
        data.cluster_multiplicity_rms = h->GetRMS();
        delete h;
      } else {
        std::cout << "run number: " << runnumber << " h_MicromegasClusterQA_cluster_multiplicity not found." << std::endl;
      }

      if(h_size)
      {
        auto h = h_size->ProjectionY( "_proj", i+1, i+1 );
        data.cluster_size_mean = h->GetMean();
        data.cluster_size_rms = h->GetRMS();
        delete h;
      } else {
        std::cout << "run number: " << runnumber << " h_MicromegasClusterQA_cluster_size not found." << std::endl;
      }

      if(h_charge)
      {
        auto h = h_charge->ProjectionY( "_proj", i+1, i+1 );
        data.cluster_charge_mean = h->GetMean();
        data.cluster_charge_rms = h->GetRMS();
        delete h;
      } else {
        std::cout << "run number: " << runnumber << " h_MicromegasClusterQA_cluster_charge not found." << std::endl;
      }
    }
  }

  std::cout << "CheckRawDataMicromegasCluster_QA - data saved" << std::endl;

  // create per detector TGraphErrors
  using tg_array_t = std::array<TGraphErrors*, 16>;
  tg_array_t tg_mult_mean;
  tg_array_t tg_size_mean;
  tg_array_t tg_charge_mean;
  for( int i = 0; i < 16; ++i )
  {
    tg_mult_mean[i] = new TGraphErrors();
    tg_size_mean[i] = new TGraphErrors();
    tg_charge_mean[i] = new TGraphErrors();

    for( auto tg:{tg_mult_mean[i], tg_size_mean[i], tg_charge_mean[i]} )
    {

      tg->SetLineColor(1);
      tg->SetLineWidth(2);
      tg->SetMarkerStyle(20);
      tg->SetMarkerColor(1);
      tg->SetMarkerSize(1);
      tg->GetXaxis()->SetTitle( "Run number" );
    }
    tg_mult_mean[i]->GetYaxis()->SetTitle("#LTN_{cluster}/trigger#GT");
    tg_size_mean[i]->GetYaxis()->SetTitle("#LTcluster size#GT (strips)");
    tg_charge_mean[i]->GetYaxis()->SetTitle("#LTcluster charge#GT (ADC)");
  }
  std::cout << "CheckRawDataMicromegasCluster_QA - tgrapherrors created" << std::endl;

  // loop over data and fill TGraphs
  std::array<int,16> ipt={};
  for( const auto& [runnumber, data_array]:data_map )
  {
    for( int i = 0; i < 16; ++i )
    {

//       if( i==0 )
//       { std::cout << "runnumber: " << runnumber << " ref: " << data_array[i].clusters_ref << std::endl; }

      if( data_array[i].entries>500 )
      {
        tg_mult_mean[i]->SetPoint(ipt[i], runnumber, data_array[i].cluster_multiplicity_mean);
        tg_size_mean[i]->SetPoint(ipt[i], runnumber, data_array[i].cluster_size_mean);
        tg_charge_mean[i]->SetPoint(ipt[i], runnumber, data_array[i].cluster_charge_mean);
        ++ipt[i];
      }
    }
  }

  for( int i = 0; i < 16; ++i )
  { std::cout << "detector: " << detector_names[i] << " entries: " << ipt[i] << std::endl; }

  std::cout << "CheckRawDataMicromegasCluster_QA - tgrapherrors filled" << std::endl;

  // make plots
  if( true )
  {
    int icv = 0;
    for( const auto& tg_array:{tg_mult_mean, tg_size_mean, tg_charge_mean} )
    {
      const auto cv_name = Form( "cv_%i", icv++ );
      auto cv = new TCanvas( cv_name, "", 1200, 800 );
      cv->Divide(4,4);
      for( int i = 0; i < 16; ++i )
      {
        cv->cd(i+1);
        gPad->SetTopMargin( 0.01 );
        gPad->SetLeftMargin( 0.1 );
        gPad->SetRightMargin( 0.01 );
        tg_array[i]->Draw("APL");

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

        Draw::PutText( 0.15, 0.85, detector_names[i].c_str(), 0.1 );
      }
      pdfDocument.Add(cv);
    }
  }


  if( true )
  {
    int icv = 0;
    for( const auto& tg_array:{tg_mult_mean, tg_size_mean, tg_charge_mean} )
    {
      const auto cv_name = Form( "cv_%i", icv++ );
      auto cv = new TCanvas( cv_name, "", 1200, 800 );
      {
        int i = 5;
        gPad->SetTopMargin( 0.01 );
        gPad->SetLeftMargin( 0.12 );
        gPad->SetRightMargin( 0.01 );
        tg_array[i]->Draw("APL");

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

        Draw::PutText( 0.15, 0.85, detector_names[i].c_str(), 0.1 );
      }
      pdfDocument.Add(cv);
    }
  }

}
