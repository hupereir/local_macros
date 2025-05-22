#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Simulations";

void DrawOccupancy()
{
  set_style( false );

  const TString rootfilename = "Rootfiles/Occupancy_sHijing_0_20fm_bkg_0_20fm-0000000007.root";
  const TString pdfFile = "Figures/Occupancy_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000007.pdf";
  PdfDocument pdfDocument( pdfFile );

  auto f = TFile::Open( rootfilename );

  const TString label[] = {
    "#Phi strips",
    "#it{Z} strips"
  };

  for( int i =0; i< 2; ++i )
  {
    auto h = static_cast<TH1*>( f->Get(Form( "h_%i", i) ) );
    h->StatOverflows(true);
    h->GetXaxis()->SetTitle( "hit count" );
    h->GetYaxis()->SetTitle( "A.U." );
    h->GetYaxis()->SetTitleOffset( 1.4 );
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );

    const auto mean = 100.*h->GetMean()/256;
    std::cout << "Occupancy - mean: " << mean << std::endl;

    {
      TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
      h->Scale( 1./h->GetEntries() );
      h->Draw("hist");
     gPad->SetLeftMargin( 0.14 );
     gPad->SetLogy(true);

      auto text = new TPaveText(0.36,0.67,0.94,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV, 50kHz pile-up");
      text->AddText(Form( "%s,  occupancy=%.1f %%", label[i].Data(), mean ) );
      // text->AddText("14/11/2023");
      text->Draw();
      pdfDocument.Add(cv);
    }

    {
      TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
      auto p = static_cast<TProfile*>( f->Get(Form("p_%i", i ) ) );
      p->GetYaxis()->SetTitle( "Occupancy (%)" );
      p->SetMaximum( 40 );
      p->SetMarkerSize(2);
      p->Draw();

      auto text = new TPaveText(0.36,0.67,0.94,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(1010);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV, 50kHz pile-up");
      text->AddText(Form( "%s,  occupancy=%.1f %%", label[i].Data(), mean ) );
      // text->AddText("14/11/2023");
      text->Draw();
      pdfDocument.Add(cv);
    }
  }

}
