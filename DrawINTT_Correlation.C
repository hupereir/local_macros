#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Preliminary";

void DrawINTT_Correlation()
{
  set_style( false );
  
  const int intt_id = 0;
  const TString pdffilename = Form( "Figures/INTT_Correlation_intt%i-all.pdf", intt_id );
  const TString pngfilename = Form( "Figures/INTT_Correlation_intt%i-all.png", intt_id );
  PdfDocument pdfDocument( pdffilename );
  
  TH2* h_correlation = nullptr;
  for( const int runnumber: {20445, 20446} )
  {
    const TString rootfilename = Form( "Rootfiles/INTT_Correlation-%08i-0000_intt%i.root", runnumber, intt_id );
    auto tfile = std::unique_ptr<TFile>( TFile::Open( rootfilename ) );
    auto h2 = static_cast<TH2*>( tfile->Get( "h_correlation" ) );
    
    if( !h_correlation )
    {
      gROOT->cd();
      h_correlation = static_cast<TH2*>( h2->Clone( "h_correlation_all" ) );
    } else {
      h_correlation->Add(h2);
    }
  }
  
  // make nice plot
  
  if( true )
  {
    TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
    h_correlation->GetXaxis()->SetTitle( "TPOT clusters" );
    h_correlation->GetYaxis()->SetTitle( Form( "INTT hits (INTT%i)", intt_id ));
    h_correlation->Draw("colz");
    h_correlation->GetYaxis()->SetTitleOffset( 1.8 );
    gPad->SetLeftMargin( 0.18);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
    
    // auto text = new TPaveText(0.5,0.16,0.86,0.34, "NDC" );
    auto text = new TPaveText(0.57,0.16,0.86,0.31, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(1010);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText("07/11/2023");
    text->Draw();
    pdfDocument.Add(cv);
    cv->SaveAs(pngfilename);
  }  
}
