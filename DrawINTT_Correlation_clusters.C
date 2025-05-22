#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Preliminary";

void DrawINTT_Correlation_clusters()
{
  set_style( false );
  const int intt_id = 0;
  const TString pdffilename = Form( "Figures/INTT_Correlation_clusters_intt%i-all.pdf", intt_id );
  const TString pngfilename = Form( "Figures/INTT_Correlation_clusters_intt%i-all.png", intt_id );
  
  TH2* h_correlation = nullptr;
  for( const int runnumber: {20445, 20446} )
  {
    const TString rootfilename = Form( "Rootfiles/INTT_Correlation_clusters-%08i-0000_intt%i.root", runnumber, intt_id );
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
    auto cv = new TCanvas( "cv", "cv", 980, 900 );
    h_correlation->GetXaxis()->SetTitle( "TPOT clusters" );
    h_correlation->GetYaxis()->SetTitle( Form("INTT clusters (INTT%i)",intt_id) );
    h_correlation->GetYaxis()->SetTitleOffset( 1.6 );
    h_correlation->Draw("colz");
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.15);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
    
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      text->DrawLatex( 0.65, 0.95, "#it{07/11/2023}" );
    }
    
    {
      auto text = new TPaveText(0.5,0.14,0.85,0.25, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->Draw();
    }
    cv->SaveAs(pdffilename);
    cv->SaveAs(pngfilename);
  }  
}
