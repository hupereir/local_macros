#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>

#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasMapping.h>

R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{
  MicromegasMapping mapping;
}

void CompareCalibrations()
{
//   PdfDocument pdfDocument( "Figures/CompareCalibrations_39495.pdf" );

//   static constexpr int nfiles = 3;
//   const TString filenames[] =
//   {
//     "Calibrations/TPOT_Pedestal-00009416-0000.root",
//     "Calibrations/TPOT_Pedestal-00039495-0000.root",
//     "Calibrations/TPOT_Pedestal-00045955-0000-from_json.root"
//   };
//
//   const TString labels[] { "run 9416", "run 39495", "run 45955 (ZS)" };
//   int colors[] = {1, 2, 4};

//   PdfDocument pdfDocument( "Figures/CompareCalibrations_45955.pdf" );
//
//   static constexpr int nfiles = 2;
//   const TString filenames[] =
//   {
//     "Calibrations/TPOT_Pedestal-00045955-0000.root",
//     "Calibrations/TPOT_Pedestal-00045955-0000-from_json.root"
//   };
//   const TString labels[] { "run 45955 (06/18/24)", "run 52122 (08/27/24)" };
//   int colors[] = {1, 2};

//   PdfDocument pdfDocument( "Figures/CompareCalibrations_52122.pdf" );
//
//   static constexpr int nfiles = 2;
//   const TString filenames[] =
//   {
//     "Calibrations/TPOT_Pedestal-00045955-0000.root",
//     "Calibrations/TPOT_Pedestal-00052122-0000.root"
//   };
//
//   const TString labels[] { "run 45955 (06/18/24)", "run 52122 (08/27/24)" };
//   int colors[] = {1, 2};

  PdfDocument pdfDocument( "Figures/CompareCalibrations_52122.pdf" );

  static constexpr int nfiles = 2;
  const TString filenames[] =
  {
    "Calibrations/TPOT_Pedestal-00039495-0000.root",
    "Calibrations/TPOT_Pedestal-00052123-0000.root"
  };

  const TString labels[] { "run 39495", "run 52123 (08/27/24)" };
  int colors[] = {1, 2};
  int symbols[] = {20, 20};

  std::vector<MicromegasCalibrationData> calibrations;
  for( const auto& filename:filenames )
  {
    MicromegasCalibrationData calibration;
    calibration.read(filename.Data());
    calibrations.push_back( calibration );
  }


  // produce plot
  auto cv = new TCanvas( "cv_tile_all", "cv_all", 1200, 1200 );
  cv->Divide( 4, 4, 0.002, 0.002 );

  const auto fee_id_list = mapping.get_fee_id_list();
  for( int ifee = 0; ifee<fee_id_list.size(); ++ifee )
  {
    cv->cd(ifee+1);

    const auto& fee_id( fee_id_list[ifee] );
    const auto hname = Form("h%i", fee_id );
    auto h = new TH1F( hname, "", 256, 0, 256 );
    h->SetStats(false);
    h->SetMinimum(0);
    h->SetMaximum(200);
    h->GetXaxis()->SetTitle("channel");
    h->GetYaxis()->SetTitle("threshold (ADC)");
    h->Draw();

    auto legend = new TLegend( 0.15, 0.19, 0.72, 0.38, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    // create TGraphs
    for( int ifile = 0; ifile < nfiles; ++ifile )
    {
      const auto& calibration( calibrations[ifile] );

      auto tge = new TGraph();
      tge->SetMarkerStyle(symbols[ifile]);
      tge->SetMarkerColor(colors[ifile]);
      for( int ich = 0; ich < 256; ++ich )
      {
        const auto pedestal = calibration.get_pedestal(fee_id, ich);
        const auto rms = calibration.get_rms(fee_id, ich);
        const auto threshold = pedestal+5.*rms;
        tge->SetPoint(ich, ich, threshold );
      }

      tge->Draw("P");
      legend->AddEntry( tge, labels[ifile], "P" );
    }

    // add detector name
    const auto detname = mapping.get_detname_sphenix(fee_id);
    Draw::PutText( 0.15, 0.8, detname.c_str(), 0.1);

    if( ifee == 0 ) legend->Draw();

  }

  pdfDocument.Add(cv);

}
