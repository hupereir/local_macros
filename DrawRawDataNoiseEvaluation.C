#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <micromegas/MicromegasMapping.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)


namespace
{
  MicromegasMapping mapping;

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

  std::vector<std::string> get_detector_names()
  {
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
    return detector_names;
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

}

//___________________________________________________________
void DrawRawDataNoiseEvaluation()
{

  gStyle->SetOptStat(0);

  const auto detector_names = get_detector_names();

  const int n_files = 2;
  const TString filenames[] = {
    "Rootfiles/RawDataNoiseEvaluation-00034567-0000.root",
      // "Rootfiles/RawDataNoiseEvaluation-00029911-0000.root",
    "Rootfiles/RawDataNoiseEvaluation-00025475-0000.root"
  };

  const int colors[] = {1, 2};
  const TString labels[] = {
    "Mar 18, Run 34567",
//    "Feb 13, Run 29904",
    "Reference, Run 25475"
  };

  auto cv = new TCanvas( "cv1", "cv1", 980, 900 );
  gPad->SetTopMargin( 0.07 );
  gPad->SetLeftMargin( 0.14);

  TH1* h = new TH1I( "dummy", "", 4096, 0, 4096 );
  h->SetMinimum(0);
  h->SetMaximum(40);
  h->GetXaxis()->SetTitle( "strip" );
  h->GetYaxis()->SetTitle( "RMS (adc)" );
  h->GetYaxis()->SetTitleOffset( 1.3 );
  h->Draw();

  auto legend = new TLegend( 0.15, 0.81, 0.62, 0.91, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  for( int i=0; i<n_files; ++i )
  {
    auto f = TFile::Open( filenames[i] );
    auto h = static_cast<TH1*>( f->Get("h_rms") );
    h->SetLineColor( colors[i] );
    h->Draw("same");
    legend->AddEntry( h, labels[i], "L" );
  }

  gPad->Update();

  draw_vertical_lines( gPad );
  draw_detector_names( gPad, detector_names, 0.02 );

  auto line = Draw::HorizontalLine( gPad, 10 );
  line->SetLineColor(2);
  line->Draw();

  cv->SaveAs("Figures/RawDataNoiseEvaluation_comparison1.pdf" );

}
