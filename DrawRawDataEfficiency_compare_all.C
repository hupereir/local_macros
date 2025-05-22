#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

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

  int get_color( int layer, int tile, int region )
  {
    std::array<int, 8> color_base = {kOrange, kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan };
    if( layer == 0 ) return color_base[tile]+region;
    else return color_base[tile]-1-region;
  }
}

//______________________________________________________
void DrawRawDataEfficiency_compare_all()
{

  set_style(false);
  save_detector_names();

  std::vector<TString> inputfiles = {
    "Rootfiles/RawDataEfficiency_new_z_2024_07_08-D400.root",
    "Rootfiles/RawDataEfficiency_new_z_2024_07_08-D400_new_calib_zs.root"
  };

//   std::vector<TString> inputfiles = {
//     "Rootfiles/RawDataEfficiency_new_phi_2024_07_08-D400.root",
//     "Rootfiles/RawDataEfficiency_new_phi_2024_07_08-D400_new_calib_zs.root"
//   };

  int symbols[] = { 24, 20 };
  TString labels[] =
  {
    "offline threshold",
    "ZS threshold"
  };


  TLegend* legends[16] = {nullptr};

  // make plot
  if( true )
  {
    auto cv = new TCanvas( "cv_tile_all", "cv_all", 1200, 1200 );
    cv->Divide( 4, 4, 0.002, 0.002 );

    // one canvas per region
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& detname = detector_names[detid];

      cv->cd( detid+1 );

      const auto hname = Form( "h_tile_%i", detid );
      auto h = new TH1I( hname, "", 100, 310, 460 );
      h->SetStats(0);
      h->SetMinimum( 0 );
      h->SetMaximum( 1. );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "efficiency" );
      h->Draw();

      legends[detid] = new TLegend( 0.15, 0.73, 0.81, 0.91, "", "NDC" );
      legends[detid]->SetFillColor(0);
      legends[detid]->SetFillStyle(0);
      legends[detid]->SetBorderSize(0);
      legends[detid]->Draw();

    }

    for( int i=0; i<inputfiles.size(); ++i )
    {
      const auto& inputfile = inputfiles[i];
      auto tfile = std::make_unique<TFile>( inputfile, "READ" );

      // one canvas per region
      for( int ilayer = 0; ilayer <2; ++ilayer )
        for( int itile = 0; itile <8; ++itile )
      {
        const int detid = itile + 8*ilayer;
        const auto& detname = detector_names[detid];
        cv->cd( detid+1 );

        auto tg = static_cast<TGraphErrors*>( tfile->Get( Form( "tge_%s", detname.c_str() ) ) );
        tg->SetMarkerStyle( symbols[i] );
        tg->Draw("P");

        legends[detid]->AddEntry( tg, Form( "%s - %s", detname.c_str(), labels[i].Data() ), "P" );
      }

    }
  }

}
