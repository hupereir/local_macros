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


//___________________________________________________________________-
void DrawRawDataClusterCharge_new()
{
  using run_map_t = std::vector<std::pair<int,int>>;
  save_detector_names();

  // field off scan (D100)
  run_map_t run_map =
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

  const TString rootfilename = "Rootfiles/RawDataClusterCharge_all_old-D100.root";
  const TString pdffilename = "Figures/RawDataClusterCharge_all_old-D100.pdf";

//   // phi scan (D400)
//   run_map_t run_map =
//   {
//     { 430, 21079 },
//     { 420, 21080 },
//     { 410, 21081 },
//     { 400, 21082 },
//     { 390, 21083 },
//     { 380, 21084 },
//     { 370, 21153 },
//     { 360, 21154 },
//     { 350, 21155 },
//     { 340, 21156 },
//     { 330, 21157 },
//     { 320, 21158 }
//   };
//
//   const TString rootfilename = "Rootfiles/RawDataClusterCharge_phi_new0-D400.root";
//   const TString pdffilename = "Figures/RawDataClusterCharge_phi_new0-D400.pdf";

//   // phi scan (D400)
//   run_map_t run_map =
//   {
//     { 430, 24341 },
//     { 420, 24342 },
//     { 410, 24343 },
//     { 400, 24347 },
//     { 390, 24348 },
//     { 380, 24349 },
//     { 370, 24350 },
//     { 360, 24351 },
//     { 350, 24352 },
//     { 340, 24353 },
//     { 330, 24354 },
//     { 320, 24355 }
//   };
//
//   const TString rootfilename = "Rootfiles/RawDataClusterCharge_phi_new-D400.root";
//   const TString pdffilename = "Figures/RawDataClusterCharge_phi_new-D400.pdf";

  std::cout << "DrawRawDataClusterCharge_new - rootfilename: " << rootfilename << std::endl;
  std::cout << "DrawRawDataClusterCharge_new - pdffilename: " << pdffilename << std::endl;

  PdfDocument pdfdocument( pdffilename );
  RootFile rootfile( rootfilename );

  // create TGraphErrors
  std::array<TGraphErrors*,16> tge_tile_array = {};
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const auto& detname = detector_names[detid];

    const auto tgname = Form( "tge_%s", detname.c_str());
    const auto tgtitle =  Form( "Cluster charge %s", detname.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("cluster charge (ADC)");
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
    const auto tgtitle =  Form( "Cluster charge %s", name.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("cluster charge (ADC)");
    tge->SetMarkerStyle( 20 );

    tge->SetMarkerColor( 1 );
    tge->SetLineColor(  1 );

    tg_mean.push_back( tge );
    rootfile.Add(tge);
  }

  // loop over run numbers
  int ipoint = 0;
  for( const auto& [hv,runnumber]:run_map )
  {
    const TString inputFile = Form( "Rootfiles/RawDataClusterCharge-%08i-0000.root", runnumber );
    std::cout << "DrawRawDataClusterCharge_new - inputFile: " << inputFile << std::endl;

    std::unique_ptr<TFile> input(TFile::Open( inputFile));

    // loop over detectors
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto hname = Form( "h2d_%i", detid );
      const auto h = static_cast<TH1*>( input->Get(hname) );
      const auto mean = h->GetMean();
      tge_tile_array[detid]->SetPoint( ipoint, hv, mean );
    }

    {
      auto h = static_cast<TH1*>( input->Get("h_phi") );
      const auto mean = h->GetMean();
      tg_mean[0]->SetPoint( ipoint, hv, mean );
    }

    {
      auto h = static_cast<TH1*>( input->Get("h_z") );
      const auto mean = h->GetMean();
      tg_mean[1]->SetPoint( ipoint, hv, mean );
    }

    ++ipoint;
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
      h->GetYaxis()->SetTitle( "cluster charge" );
      h->SetMinimum(100);
      h->SetMaximum(2000);
      h->Draw();

      tge_tile_array[detid]->Draw("P");
      Draw::PutText( 0.15, 0.8, detname.c_str(), 0.1);
      gPad->SetLogy( true );
    }

    pdfdocument.Add( cv );

  }

}
