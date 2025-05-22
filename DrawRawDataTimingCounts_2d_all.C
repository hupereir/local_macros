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

}

//_________________________________________________
void DrawRawDataTimingCounts_2d_all()
{

  gStyle->SetOptStat(0);

  set_style( false );

  using run_pair_t = std::pair<int,int>;
  using run_map_t = std::vector<run_pair_t>;

  run_map_t run_map =
  {
    {50, 20464},
    {100, 20465},
    {150, 20466},
    {200, 20467},
    {250, 20468},
    {300, 20469},
    {350, 20470},
    {400, 20471}
  };

  save_detector_names();

  // make plot
  const TString pdffilename = "Figures/RawDataTimingCounts_2d_drift.pdf";
  PdfDocument pdfDocument( pdffilename );

  // loop over runs
  int ipoint = 0;
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
//   int ilayer = 0;
//   int itile = 0;
  {

    // create canvas
    const int detid = itile + 8*ilayer;
    auto cv = new TCanvas( Form( "cv_%i", detid ), "", 1200, 600 );
    cv->Divide( 4, 2, 0.001, 0.001 );

    int icv = 0;
    for( const auto& [hv,runnumber]:run_map )
    {

      const TString inputfilename = Form( "Rootfiles/RawDataTimingCounts_2d-%08i-0000.root", runnumber );
      auto tfile = std::unique_ptr<TFile>( TFile::Open( inputfilename, "READ") );

      const auto hname = Form( "h_%i_%i", ilayer, itile );
      const auto h = static_cast<TH2*>( tfile->Get( hname ) );
      h->GetYaxis()->SetTitleOffset(1.5);
      h->SetTitle( "" );

      cv->cd( ++icv );

      // draw
      h->SetMaximum(std::min<double>(100, h->GetMaximum() ));
      h->DrawCopy( "colz" );
      gPad->SetRightMargin( 0.2 );
      gPad->SetLeftMargin( 0.15 );

      auto text = new TPaveText(0.16,0.76,0.96,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "#it{#bf{sPHENIX}}" );
      text->AddText( Form( "%s, V_{drift}=-%iV", detector_names[detid].c_str(), hv ) );
      text->Draw();
    }

    // cv->SaveAs(Form("Figures/RawDataTimingCounts_2d_drift_%s.png", detector_names[detid].c_str() ) );
    cv->SaveAs(Form("Figures/RawDataTimingCounts_2d_drift_%s.pdf", detector_names[detid].c_str() ) );
    pdfDocument.Add( cv );
  }
}
