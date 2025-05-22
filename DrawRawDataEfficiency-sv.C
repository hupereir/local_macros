#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

const std::string status = "Internal";

void DrawRawDataEfficiency()
{
  const TString inputfile = "Rootfiles/RawDataEfficiency_phi-D400_4sigma_charge_size_cut.root";
  auto tfile = std::make_unique<TFile>( inputfile, "READ" );
  
  const TString detector = "SCIP"; 
  
  auto tg = static_cast<TGraphErrors*>( tfile->Get( Form( "tge_%s", detector.Data() ) ) );
  tg->SetTitle( "" );
  
  
  // create plot
  TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
  tg->SetMarkerSize( 2 );
  tg->Draw("AP" );
  gPad->SetTopMargin( 0.07 );
  gPad->SetLeftMargin( 0.14);

  {
    auto text = new TLatex;
    text->SetNDC( true );
    text->SetTextSize( 0.045 );
    text->DrawLatex( 0.77, 0.95, "#it{07/31/2023}" );
  }

  {
    auto text = new TPaveText(0.17,0.73,0.53,0.88, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText(Form( "Detector: %s", detector.Data() ) );
    text->Draw();
  }
}
