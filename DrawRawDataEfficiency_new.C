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

//____________________________________________________________________
void DrawRawDataEfficiency_new()
{

  save_detector_names();

//   const TString inputfile = "Rootfiles/RawDataEfficiency_new_phi_2024_08_09-D400.root";
//   const TString tag = "_phi_2024_08_09-D400";
//   int ilayer = 0;

  const TString inputfile = "Rootfiles/RawDataEfficiency_new_z_2024_08_09-D400.root";
  const TString tag = "_z_2024_08_09-D400";
  int ilayer = 1;

//   const TString inputfile = "Rootfiles/RawDataEfficiency_new_phi_2024_06_27-D400.root";
//   const TString tag = "_phi_2024_06_27-D400";
//   int ilayer = 0;

//   const TString inputfile = "Rootfiles/RawDataEfficiency_new_z_2024_07_08-D400.root";
//   const TString tag = "_z_2024_07_08-D400";
//   int ilayer = 1;

  const TString pdffile = Form( "Figures/RawDataEfficiency_compare%s.pdf", tag.Data() );

  auto tfile = TFile::Open( inputfile, "READ" );

  // create plot
  TCanvas* cv0 = new TCanvas( "cv0", "cv", 980, 900 );
  {
    auto h = new TH1I( "h", "", 100, 310, 460 );
    h->SetStats(false);
    h->SetMinimum(0);
    h->SetMaximum(1);
    h->GetXaxis()->SetTitle( "resist HV (V)" );
    h->GetYaxis()->SetTitle( "efficiency" );
    h->Draw();
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    gPad->Update();
  }

  auto legend = new TLegend( 0.2, 0.5, 0.5, 0.8, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  // loop over modules
  // for( int itile = 0; itile < 8; ++itile )
  int itile = 3;
  {
    const auto detector_name = detector_names[itile+8*ilayer];
    auto color = get_color(ilayer, itile, 0 );
    auto tg = static_cast<TGraphErrors*>( tfile->Get( Form( "tge_%s", detector_name.c_str() ) ) );
    tg->SetMarkerSize(2);
    tg->Draw("P");

    legend->AddEntry(tg, detector_name.c_str(), "PL" );
  }

  legend->Draw();

  cv0->SaveAs(pdffile);

}
