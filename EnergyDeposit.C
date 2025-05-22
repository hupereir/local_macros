#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//________________________________________________
TString EnergyDeposit()
{
  const TString inputFile =  "DST/CONDOR_SimEvaluation_hijing/dst_simeval_sHijing_0_20fm-0000000007-00*.root";
  const TString pdfFile = "Figures/EnergyDeposit_sHijing_0_20fm-0000000007.pdf";

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const double edep_max = 20;
  const TString var( "_edep*1e6" );
  // const TCut cut = Form( "_edep*1e6<%f && _delta_r > 0.2", edep_max );
  const TCut cut = "_delta_r > 0.2";

  auto h = new TH1F( "edep", "", 100, 0, edep_max );
  h->StatOverflows(true);
  h->GetXaxis()->SetTitle( "E_{deposited} (keV)" );
  h->GetYaxis()->SetTitle( "A.U." );
  h->SetFillStyle(1001);
  h->SetFillColor(kYellow );

  Utils::TreeToHisto( tree, h->GetName(), var, cut, false );

  const double mean = h->GetMean();
  const double mpv = h->GetBinCenter( h->GetMaximumBin() );

  std::cout << "EnergyDeposit - mean: " << mean << std::endl;
  std::cout << "EnergyDeposit - mpv: " << mpv << std::endl;

  // create canvas and divide
  auto cv = new TCanvas( "cv", "cv", 980, 900 );
  gPad->SetTopMargin( 0.05 );
  gPad->SetRightMargin( 0.05);

  h->Draw();
  gPad->SetLogy( true );

  {
    auto text = new TPaveText(0.57,0.78,0.93,0.90, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText("HIJING - Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText(Form("mean: %.2f keV", mean ) );
    text->AddText(Form("MPV: %.2f keV", mpv ) );
    text->Draw();
  }

  pdfDocument.Add(cv);
  return pdfFile;

}
