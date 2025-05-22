#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Preliminary";

void DrawTPOT_Correlation()
{
  set_style( false );
  const TString pdffilename = "Figures/TPOT_Correlation-all.pdf";
  const TString pngfilename = "Figures/TPOT_Correlation-all.png";
  
  TH2* h_correlation = nullptr;
  for( const int runnumber: {20445, 20446} )
  {
    const TString rootfilename = Form( "Rootfiles/TPOT_Correlation-%08i-0000.root", runnumber );
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
    h_correlation->Draw("colz");
    h_correlation->GetXaxis()->SetTitle( "TPOT clusters (#it{#phi} views)" );
    h_correlation->GetYaxis()->SetTitle( "TPOT clusters (#it{Z} views)" );
    h_correlation->GetYaxis()->SetTitleOffset( 1.4 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
    
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      text->DrawLatex( 0.65, 0.95, "#it{07/11/2023}" );
    }
    
    {
      auto text = new TPaveText(0.51,0.15,0.87,0.26, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->Draw();
    }
    
    cv->SaveAs(pdffilename);
    cv->SaveAs(pngfilename);
  }  
}
