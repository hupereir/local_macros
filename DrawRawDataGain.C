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

  TF1* fit_tgraph( TGraph* tg, double xmin = 390, double xmax = 420 )
  {

    const auto fname = Form( "%s_fit", tg->GetName() );
    auto f = new TF1( fname, "expo", 320, 450 );
    tg->Fit( f, "0Q", "", xmin, xmax );
    f->SetLineColor(2);
    f->Draw("same");
    return f;
  }

}

//___________________________________________________________________-
void DrawRawDataGain()
{
  using run_map_t = std::vector<std::pair<int,int>>;
  save_detector_names();

  // conversion factors
  const double electron_charge = 1.602e-19 * 1e15; // fC
  const double fee_gain_phi = 16; // mv/fC
  const double fee_gain_z = 13.5; // mV/fc
  const double adc_conv = 2.15; // mV/ADC
  const double n_primary = 1.26e3 / 26.3; // electrons

  const double scale_phi = adc_conv/(fee_gain_phi*electron_charge*n_primary);
  const double scale_z = adc_conv/(fee_gain_z*electron_charge*n_primary);

  #if false
  // field off scan (D100)
  const run_map_t run_map_phi =
  {
    { 430, 9439 },
    { 420, 9440 },
    { 410, 9441 },
    { 400, 9443 },
    { 390, 9444 },
    { 380, 9446 },
    { 370, 9447 },
    { 360, 9448 },
    { 350, 9449 },
    { 340, 9450 },
    { 330, 9451 },
    { 320, 9452 }
  };

  const run_map_t run_map_z = run_map_phi;

  const TString rootfilename = "Rootfiles/RawDataGain_all_old-D100.root";
  const TString pdffilename = "Figures/RawDataGain_all_old-D100.pdf";

  #endif

  #if false
  // phi scan (D400)
  const run_map_t run_map_phi =
  {
    { 430, 21079 },
    { 420, 21080 },
    { 410, 21081 },
    { 400, 21082 },
    { 390, 21083 },
    { 380, 21084 },
    { 370, 21153 },
    { 360, 21154 },
    { 350, 21155 },
    { 340, 21156 },
    { 330, 21157 },
    { 320, 21158 }
  };

  // z scan (D400)
  const run_map_t run_map_z =
  {
    { 430, 21159 },
    { 420, 21160 },
    { 410, 21161 },
    { 400, 21162 },
    { 390, 21163 },
    { 380, 21164 },
    { 370, 21166 },
    { 360, 21168 },
    { 350, 21169 },
    { 340, 21170 },
    { 330, 21171 },
    { 320, 21172 }
  };

  const TString rootfilename = "Rootfiles/RawDataGain_new0-D400.root";
  const TString pdffilename = "Figures/RawDataGain_new0-D400.pdf";
  #endif

  #if true
  // phi scan (D400)
  const run_map_t run_map_phi =
  {
    { 430, 24341 },
    { 420, 24342 },
    { 410, 24343 },
    { 400, 24347 },
    { 390, 24348 },
    { 380, 24349 },
    { 370, 24350 },
    { 360, 24351 },
    { 350, 24352 },
    { 340, 24353 },
    { 330, 24354 },
    { 320, 24355 }
  };

  // z scan (D400)
  const run_map_t run_map_z =
  {
    {430, 24358},
    {420, 24359},
    {410, 24360},
    {400, 24365},
    {390, 24366},
    {380, 24367},
    {370, 24368},
    {360, 24369},
    {350, 24370},
    {340, 24371},
    {330, 24372},
    {320, 24373}
  };


  const TString rootfilename = "Rootfiles/RawDataGain_new-D400.root";
  const TString pdffilename = "Figures/RawDataGain_new-D400.pdf";
  #endif

  std::cout << "DrawRawDataGain_new - rootfilename: " << rootfilename << std::endl;
  std::cout << "DrawRawDataGain_new - pdffilename: " << pdffilename << std::endl;

  PdfDocument pdfdocument( pdffilename );
  RootFile rootfile( rootfilename );

  static const double nominal_hv = 400;
  static const double max_hv = 450;
  std::array<double,16> nominal_gain = {};
  double nominal_gain_phi = 0;
  double nominal_gain_z = 0;


  // create TGraphErrors
  std::array<TGraphErrors*,16> tge_tile_array = {};
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const auto& detname = detector_names[detid];

    const auto tgname = Form( "tge_%s", detname.c_str());
    const auto tgtitle =  Form( "Gain %s", detname.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("Gain");
    tge->SetMarkerStyle( 20 );

    tge->SetMarkerColor( get_color( ilayer, itile, 0 ) );
    tge->SetLineColor( get_color( ilayer, itile, 0 ) );

    tge_tile_array[detid]=tge;
    rootfile.Add(tge);
  }

  std::vector<TGraphErrors*> tg_mean;
  for( const std::string& name:{ "phi", "z"} )
  {
    const auto tgname = Form( "tge_%s", name.c_str());
    const auto tgtitle =  Form( "Gain %s", name.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("Gain");
    tge->SetMarkerStyle( 20 );
    tge->SetMarkerSize( 2 );


    tge->SetMarkerColor( 1 );
    tge->SetLineColor(  1 );

    tg_mean.push_back( tge );
    rootfile.Add(tge);
  }

  // loop over run numbers
  for( int ilayer = 0; ilayer < 2; ++ilayer )
  {

    const auto& run_map = (ilayer == 0) ? run_map_phi:run_map_z;

    int ipoint = 0;
    for( const auto& [hv,runnumber]:run_map )
    {
      const TString inputFile = Form( "Rootfiles/RawDataClusterCharge-%08i-0000.root", runnumber );
      std::cout << "DrawRawDataClusterCharge_new - inputFile: " << inputFile << std::endl;

      std::unique_ptr<TFile> input(TFile::Open( inputFile));

      // loop over detectors
      for( int itile = 0; itile < 8; ++itile )
      {

        const bool is_phi = (ilayer == 0);
        const double scale = is_phi ? scale_phi:scale_z;

        const int detid = itile + 8*ilayer;
        const auto hname = Form( "h2d_%i", detid );
        const auto h = static_cast<TH1*>( input->Get(hname) );
        const auto mean = h->GetMean();
        tge_tile_array[detid]->SetPoint( ipoint, hv, mean*scale );
        if( hv == nominal_hv ) nominal_gain[detid] = mean*scale;
      }

      if( ilayer == 0 )
      {

        auto h = static_cast<TH1*>( input->Get("h_phi") );
        const auto mean = h->GetMean();
        tg_mean[0]->SetPoint( ipoint, hv, mean*scale_phi );
        if( hv == nominal_hv ) nominal_gain_phi = mean*scale_phi;

      } else {

        auto h = static_cast<TH1*>( input->Get("h_z") );
        const auto mean = h->GetMean();
        tg_mean[1]->SetPoint( ipoint, hv, mean*scale_z );
        if( hv == nominal_hv ) nominal_gain_z = mean*scale_z;

      }
      ++ipoint;
    }
  }

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
      auto h = new TH1I( hname, "", 100, 310, 440 );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "detector gain" );
      h->SetMinimum(1e3);
      h->SetMaximum(1e5);
      h->Draw();

      tge_tile_array[detid]->Draw("P");
      {
        auto text = new TPaveText(0.2,0.65,0.9,0.84, "NDC" );
        text->SetFillColor(0);
        text->SetFillStyle(0);
        text->SetBorderSize(0);
        text->SetTextAlign(11);
        text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
        text->AddText(Form( "%s", detname.c_str()) );
        text->AddText(Form( "Gain (%.0f V) = %.0f", nominal_hv, nominal_gain[detid] ) );
        text->Draw();
      }
      gPad->SetLogy( true );
    }

    pdfdocument.Add( cv );

  }

  if( true )
  {
    auto cv = new TCanvas( "cv_phi", "cv_all", 1200, 1200 );
    const auto hname = "h_phi";
    auto h = new TH1I( hname, "", 100, 310, 440 );
    h->GetXaxis()->SetTitle( "Resist HV (V)" );
    h->GetYaxis()->SetTitle( "detector gain" );
    h->SetMinimum(1e3);
    h->SetMaximum(1e5);
    h->Draw();

    tg_mean[0]->Draw("P");

    auto f = fit_tgraph( tg_mean[0] );
    auto max_gain = f->Eval( max_hv );

    {
      auto text = new TPaveText(0.2,0.65,0.9,0.84, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText("#Phi strips" );
      text->AddText(Form( "Gain (%.0fV) = %.0f", nominal_hv, nominal_gain_phi ) );
      text->AddText(Form( "Gain (%.0fV) = %.0f", max_hv, max_gain ) );
      text->Draw();
    }

    gPad->SetLogy( true );
    pdfdocument.Add( cv );

  }

  if( true )
  {
    auto cv = new TCanvas( "cv_z", "cv_all", 1200, 1200 );
    const auto hname = "h_z";
    auto h = new TH1I( hname, "", 100, 310, 440 );
    h->GetXaxis()->SetTitle( "Resist HV (V)" );
    h->GetYaxis()->SetTitle( "detector gain" );
    h->SetMinimum(1e3);
    h->SetMaximum(1e5);
    h->Draw();

    tg_mean[1]->Draw("P");

    auto f = fit_tgraph( tg_mean[1] );
    auto max_gain = f->Eval( max_hv );

    {
      auto text = new TPaveText(0.2,0.65,0.9,0.84, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText("z strips" );
      text->AddText(Form( "Gain (%.0fV) = %.0f", nominal_hv, nominal_gain_z ) );
      text->AddText(Form( "Gain (%.0fV) = %.0f", max_hv, max_gain ) );
      text->Draw();
    }

    gPad->SetLogy( true );
    pdfdocument.Add( cv );

  }
}
