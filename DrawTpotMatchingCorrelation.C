#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

const std::string status = "Internal";

void DrawTpotMatchingCorrelation()
{
  set_style( false );
  const TString pdffilename = "Figures/TpotMatchingCorrelation.pdf";
  const TString pngfilename = "Figures/TpotMatchingCorrelation.png";
  PdfDocument pdfDocument( pdffilename );

  const TString rootfilename = "Rootfiles/TpotMatchingCorrelation_CombinedDataReconstruction_zf_corrected-0_00052077-00052078.root";

  auto tfile = TFile::Open( rootfilename );

  auto h_correlation = static_cast<TH2*>(tfile->Get("h_phi_all"));

  TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
  h_correlation->GetXaxis()->SetTitle( "TPOT cluster position (cm)" );
  h_correlation->GetYaxis()->SetTitle( "TPC cluster position (cm)" );
  h_correlation->Draw("col");
  // h_correlation->GetYaxis()->SetTitleOffset( 1.8 );
  gPad->SetLeftMargin(0.14);
  gPad->SetTopMargin(0.118);
  // gPad->SetLogz(true);

//   auto text = new TPaveText(0.5,0.16,0.86,0.34, "NDC" );
//   text->SetFillColor(0);
//   text->SetFillStyle(0);
//   // text->SetFillStyle(1010);
//   text->SetBorderSize(0);
//   text->SetTextAlign(11);
//   text->AddText(Form("#it{#bf{sPHENIX}} %s", status.c_str()));
//   text->AddText("p+p #sqrt{s_{NN}} = 200 GeV");
//   text->AddText("08/26/2024");
//   text->Draw();

  auto text = new TPaveText(0.1,0.8,0.98,1., "NDC" );
  text->SetFillColor(0);
  text->SetFillStyle(0);
  // text->SetFillStyle(1010);
  text->SetBorderSize(0);
  text->SetTextAlign(11);
  text->AddText(Form("#it{#bf{sPHENIX}} %s, p+p #sqrt{s_{NN}} = 200 GeV, 08/26/2024", status.c_str()));
  // text->AddText("");
  text->Draw();

  pdfDocument.Add(cv);
  cv->SaveAs(pngfilename);
}
