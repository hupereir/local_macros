#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


TGraph* DrawEfficiency_BNL( int layer )
{

  if( layer == 0 )
  {
    const int n_points  = 12;
    double x[12] = { 320, 329.92, 339.84, 349.95, 359.87, 369.79, 379.96, 389.82, 399.74, 409.85, 419.97, 429.95 };
    double y[12] = { 0.16771, 0.21759, 0.27241, 0.3289, 0.43378, 0.54095, 0.6635, 0.73649, 0.79793, 0.84616, 0.84376, 0.88052 };
    double y_err_syst_plus[12] = { 0.069358, 0.074306, 0.079969, 0.082144, 0.081377, 0.083661, 0.067325, 0.067066, 0.064322, 0.066155, 0.062325, 0.067039 };
    double y_err_syst_minus[12] = { 0.06444, 0.07764, 0.079984, 0.082225, 0.081145, 0.081687, 0.070995, 0.069068, 0.069724, 0.067892, 0.065038, 0.06912 };

    auto tg = new TGraph( n_points, x, y );
    tg->SetMarkerStyle(20);
    tg->SetMarkerColor(kMagenta);
    tg->SetLineColor(kMagenta);
    tg->SetMarkerSize(2);
    gPad->Update();

    Draw::SetBoxFillStyle(1011);
    const int merge = Color( kMagenta ).Merge( kWhite, 0.2 );
    for( int i = 0; i < n_points; ++i )
    { Draw::DrawBox( x[i], y[i], 2, 0.5*(y_err_syst_plus[i]+y_err_syst_minus[i]), merge ); }

    tg->Draw("P");

    return tg;
  }

  if( layer == 1 )
  {
    const int n_points  = 12;
    double x[12] = { 320, 329.8, 339.84, 349.87, 359.91, 369.95, 379.87, 390.03, 399.95, 410.11, 420.02, 430.06 };
    double y[12] = { 0.20135, 0.25202, 0.31887, 0.38787, 0.47844, 0.569, 0.69299, 0.7372, 0.80189, 0.81698, 0.83423, 0.84609 };
    double y_err_syst_plus[12] = { 0.079784, 0.078706, 0.079717, 0.078706, 0.079784, 0.072237, 0.070081, 0.072237, 0.070081, 0.073315, 0.071159, 0.071159 };
    double y_err_syst_minus[12] = { 0.078706, 0.081941, 0.079717, 0.077628, 0.075472, 0.065768, 0.066846, 0.073315, 0.066846, 0.069003, 0.066846, 0.074394 };

    auto tg = new TGraph( n_points, x, y );
    tg->SetMarkerStyle(20);
    tg->SetMarkerColor(kMagenta-1);
    tg->SetLineColor(kMagenta-1);
    tg->SetMarkerSize(2);

    gPad->Update();

    Draw::SetBoxFillStyle(1011);
    const int merge = Color( kMagenta-1 ).Merge( kWhite, 0.2 );
    for( int i = 0; i < n_points; ++i )
    {
      Draw::DrawBox( x[i], y[i], 2, 0.5*(y_err_syst_plus[i]+y_err_syst_minus[i]), merge );
    }

    tg->Draw("P");
    return tg;
  }

  return nullptr;

}

void DrawRawDataEfficiency_compare()
{

//   std::vector<TString> inputfiles = {
//     "Rootfiles/RawDataEfficiency_new_phi_2024_06_27-D400.root",
//       // "Rootfiles/RawDataEfficiency_new_phi_2024_07_08-D400.root",
//     "Rootfiles/RawDataEfficiency_new_phi_2024_08_09-D400.root"
//   };
//
//   const std::vector<TString> detectors =
//   {
//     "NCOP",
//     "NCOP"
//   };
//
//   int colors[] = { 1, 2 };
//   TString labels[] =
//   {
//     "NCOP - w/ ZS",
//     "NCOP - W/o ZS"
//   };

  std::vector<TString> inputfiles = {
    // "Rootfiles/RawDataEfficiency_new_phi_2024_06_27-D400.root",
    "Rootfiles/RawDataEfficiency_new_z_2024_07_08-D400.root",
    "Rootfiles/RawDataEfficiency_new_z_2024_08_09-D400.root"
  };

  const std::vector<TString> detectors =
  {
    "NCOZ",
    "NCOZ"
  };

  int colors[] = { 1, 2 };
  TString labels[] =
  {
    "NCOZ - w/ ZS",
    "NCOZ - W/o ZS"
  };

  // create plot
  TCanvas* cv0 = new TCanvas( "cv0", "cv", 980, 900 );
  {
    auto h = new TH1I( "h", "", 100, 310, 460 );
    h->SetStats(0);
    h->SetMinimum(0);
    h->SetMaximum(1);
    h->GetXaxis()->SetTitle( "resist HV (V)" );
    h->GetYaxis()->SetTitle( "efficiency" );
    h->Draw();
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    gPad->Update();
  }


  auto legend = new TLegend( 0.15, 0.76, 0.45, 0.91, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

//   {
//     auto tg = DrawEfficiency_BNL(0);
//     legend->AddEntry( tg, "NCOP - Au+Au data", "PL" );
//   }

  {
    auto tg = DrawEfficiency_BNL(1);
    legend->AddEntry( tg, "NCOZ - Au+Au data", "PL" );
  }

  for(size_t i=0; i < inputfiles.size(); ++i)
  {
    const auto& inputfile = inputfiles[i];
    const auto& detector = detectors[i];

    std::cout << "DrawRawDataEfficiency - file: " << inputfile << std::endl;
    auto tfile = std::make_unique<TFile>( inputfile, "READ" );
    auto tg = static_cast<TGraphErrors*>( tfile->Get( Form( "tge_%s", detector.Data() ) ) );

    tg->SetMarkerSize(2);
    tg->SetMarkerColor( colors[i] );
    tg->SetLineColor( colors[i] );
    tg->Draw("P");

    legend->AddEntry( tg, labels[i], "PL" );

  }

  legend->Draw();
}
