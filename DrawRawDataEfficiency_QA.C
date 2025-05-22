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
void DrawRawDataEfficiency_QA()
{
  set_style(false);
  save_detector_names();

  const auto pdfname = "Figures/RawDataEfficiency_QA_compare.pdf";
  PdfDocument pdfDocument(pdfname);

  const std::vector<TString> filenames =
  {
    "Rootfiles/RawDataEfficiency_QA.root",
    "Rootfiles/RawDataEfficiency_QA-new.root"
  };

  const std::vector<TString> labels =
  {
    "SingleMicromegasPoolInput",
    "SingleMicromegasPoolInput_v2"
  };

  const std::vector<int> colors = {1,2};

  std::vector<TFile*> tfiles;
  for( const auto& filename:filenames )
  { tfiles.push_back( TFile::Open(filename) ); }

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

      auto legend = new TLegend( 0.06, 0.18, 0.53, 0.33, "", "NDC" );
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);

      bool first = true;
      const auto tgname = Form( "tg_eff_%i", i );
      for( size_t j=0; j<tfiles.size(); ++j )
      {
        auto tg = static_cast<TGraphErrors*>(tfiles[j]->Get(tgname) );
        tg->SetMarkerColor(colors[j]);
        tg->GetXaxis()->SetRange();
        if( first ) tg->Draw("AP");
        else tg->Draw("P");

        legend->AddEntry( tg, labels[j], "P" );
        first = false;
      }

      legend->Draw();

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

}


