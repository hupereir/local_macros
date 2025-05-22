
#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

const std::string status = "Internal";

namespace
{
  
  TLatex* put_text( Double_t x, Double_t y, TString value, Double_t fontSize )
  {
    auto text = new TLatex;
    text->SetTextSize( fontSize );
    text->DrawLatex( x, y, value );
    return text;
  }

  //______________________________________________________
  void draw_vertical_lines( TVirtualPad* pad )
  {
    pad->Update();
    for( int i = 1; i < 16; ++i )
    { Draw::VerticalLine( pad, 256*i )->Draw(); }
  }

  //______________________________________________________
  void draw_detector_names( TVirtualPad* pad, std::vector<std::string> detector_names, Double_t fontSize = 0.04 )
  {
    pad->Update();
    const Double_t xMin = pad->GetUxmin();
    const Double_t xMax = pad->GetUxmax();
    const Double_t yMax = pad->GetUymax();
    const auto dx = (xMax-xMin)/16;

    for( int i = 0; i < 16; ++i )
    {
      auto x = xMin + 0.1*dx + i*dx;
      put_text( x, 0.8*yMax, detector_names[i].c_str(), fontSize )->Draw();
    }
  }
  
  //______________________________________________________
  void draw_information( int runNumber )
  {

    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      text->DrawLatex( 0.77, 0.95, "#it{09/14/2023}" );
    }
    
    {
      auto text = new TPaveText(0.17,0.73,0.53,0.88, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      text->AddText(Form( "Run %i (Cosmics)", runNumber));
      text->AddText("TPOT");
      text->Draw();
    }
  }

}

//_____________________________________________________________________________
void RawDataEvaluation(int runNumber = 25475)
{
  MicromegasMapping mapping;

  const TString inputFile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataEvaluation-%08i-0000.pdf", runNumber );
  const TString pngFile = Form( "Figures/RawDataEvaluation-%08i-0000.png", runNumber );
    
  std::cout << "RawDataEvaluation - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataEvaluation - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  
  // get 2d histogram
  // auto h_adc_channel = new TH2I( "h_adc_channel", "", 4096, 0, 4096, 1200, 0, 1200 );
  auto h_adc_channel = new TH2I( "h_adc_channel", "", 4096, 0, 4096, 500, 0, 500 );
  // const TString channel_var = "samples.channel+256*(samples.tile+8*(samples.layer-55))";
  const TString channel_var = "samples.strip+256*(samples.tile+8*(samples.layer-55))";
  const TString var = Form( "samples.adc:%s", channel_var.Data() );
  
  TCut packet_cut;
  TCut sample_cut = "samples.sample >= 5 && samples.sample <= 15";
  
  // TCut packet_cut = "samples.packet_id == 5001";
  
  Utils::max_entries = 2000;
  Utils::TreeToHisto( tree, h_adc_channel->GetName(), var, packet_cut && sample_cut, false );
  h_adc_channel->GetXaxis()->SetTitle( "strip" );
  h_adc_channel->GetYaxis()->SetTitle( "adc" );
  h_adc_channel->GetYaxis()->SetTitleOffset( 1.3 );
    
  // get detector names that match tile/layer
  std::vector<std::string> detector_names;
  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
  {
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
    const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
    detector_names.push_back( std::move( name ) );
  }
    
  // get mean and rms histogram
  auto h_pedestal = new TH1I( "h_pedestal", "", 4096, 0, 4096 );
  h_pedestal->GetXaxis()->SetTitle( "strip" );
  h_pedestal->GetYaxis()->SetTitle( "adc" );
  h_pedestal->GetYaxis()->SetTitleOffset( 1.3 );

  auto h_rms = new TH1I( "h_rms", "", 4096, 0, 4096 );
  h_rms->GetXaxis()->SetTitle( "strip" );
  h_rms->GetYaxis()->SetTitle( "adc" );
  h_rms->GetYaxis()->SetTitleOffset( 1.3 );

  TProfile* profile = nullptr;
  
  // ranged pedestals
  bool do_range = false;  
  if( do_range )
  {
    const double range_min = 0;
    const double range_max = 200;
    const auto bin_min = h_adc_channel->GetYaxis()->FindBin( range_min );
    const auto bin_max = h_adc_channel->GetYaxis()->FindBin( range_max );    
    profile = h_adc_channel->ProfileX("h_adc_channel_profx", bin_min, bin_max, "s" );
  } else {
    profile = h_adc_channel->ProfileX("h_adc_channel_profx", 1, -1, "s" );
  }
  
  for( int i = 0; i < h_adc_channel->GetNbinsX(); ++i )
  {
    h_pedestal->SetBinContent( i+1, profile->GetBinContent(i+1) );
    h_rms->SetBinContent(i+1, profile->GetBinError(i+1) );
  }
      
  {
    auto cv = new TCanvas( "cv0", "cv0", 980, 900 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);

    h_adc_channel->Draw("col");
//     draw_vertical_lines( gPad );
//     draw_detector_names( gPad, detector_names, 0.02 );
    
    draw_information( runNumber );
    
    pdfDocument.Add(cv);    
    cv->SaveAs( pngFile );
  }
  
  {
    auto cv = new TCanvas( "cv1", "cv1", 980, 900 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    
    h_pedestal->SetLineColor(1);
    h_pedestal->Draw();
   
    draw_information( runNumber );

//     draw_vertical_lines( gPad );
//     draw_detector_names( gPad, detector_names, 0.02 );
    pdfDocument.Add(cv);    

  }
  
  {
    auto cv = new TCanvas( "cv1", "cv1", 980, 900 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    
    h_rms->SetLineColor(1);
    h_rms->SetMinimum(0);
    h_rms->SetMaximum(40);
    h_rms->Draw();

    draw_information( runNumber );

//     draw_vertical_lines( gPad );    
//     draw_detector_names( gPad, detector_names, 0.02 );
    pdfDocument.Add(cv);
  }

}
