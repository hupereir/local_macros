#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Preliminary";

//_____________________________________________________________________________
TString DrawRawDataTimingCounts( int runNumber = 20126 )
{
  set_style( false );
  const TString rootfilename = Form( "Rootfiles/RawDataTimingCounts_masked-%08i-0000.root", runNumber );
  const TString pdffilename = Form( "Figures/RawDataTimingCounts-%08i-0000.pdf", runNumber );

  auto tfile = TFile::Open( rootfilename, "READ");

  int ilayer = 0;
  int itile = 0;
  int iregion = 0;

  const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
  const auto h = static_cast<TH2*>( tfile->Get( hname ) );
  h->GetXaxis()->SetTitle( "sample [50 ns]" );
  h->GetYaxis()->SetTitle( "counts" );
  h->GetYaxis()->SetTitleOffset( 2.0 );
  h->SetTitle("");

  // create canvas and divide
  auto cv = new TCanvas( "cv", "cv", 900, 900 );

  gPad->SetTopMargin( 0.07 );
  gPad->SetLeftMargin( 0.2);
  gPad->SetRightMargin( 0.05);

  h->GetXaxis()->SetRangeUser(0,100);
  h->SetFillStyle(1001);
  h->SetFillColor(kYellow );
  h->Draw("h");
  cv->Update();
  Draw::VerticalLine(gPad, 20 )->Draw();
  Draw::VerticalLine(gPad, 40 )->Draw();

  if( true )
  {
    auto text = new TLatex;
    text->SetNDC( true );
    text->SetTextSize( 0.045 );
    text->DrawLatex( 0.75, 0.95, "#it{08/30/2023}" );
  }

  {
    auto text = new TPaveText(0.57,0.78,0.93,0.90, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    // text->AddText( "#it{#bf{sPHENIX}}" );
    text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText("TPOT SCOP (#phi strips)");
    text->Draw();
  }
  cv->SaveAs(pdffilename);

  return pdffilename;

}
