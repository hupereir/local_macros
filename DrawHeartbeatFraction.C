#include <RootUtil/PdfDocument.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

// get a copy of TH1 histogram
TH1F* get_copy( TH1* h )
{
  auto copy = new TH1F(Form( "%s_norm", h->GetName() ), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  copy->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  copy->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

  for( int i = 0; i < h->GetNbinsX(); ++i )
  {
    copy->GetXaxis()->SetBinLabel(i+1, h->GetXaxis()->GetBinLabel(i+1));
    copy->SetBinContent( i+1, h->GetBinContent(i+1) );
  }

  return copy;
}


void DrawHeartbeatFraction( int runnumber=60114)
{
  const TString filename = Form("Rootfiles/Fun4All_ReadCombinedData_hp_QA-%08i-0000.root", runnumber);
  std::cout << "DrawHeartbeatFraction - filename: " << filename << std::endl;

    // create plots
  PdfDocument pdfDocument(Form("Figures/HeartbeatFraction-%08i-0000.pdf", runnumber));

  auto tfile = TFile::Open(filename);

  if( true )
  {
    auto h = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_heartbeat_stat" ) );
    auto copy = get_copy(h);

    // normalize
    auto norm = h->GetBinContent(1);

    std::cout << "DrawPackets - heartbeats: " << norm << std::endl;

    for( int i = 0; i < h->GetNbinsX(); ++i )
    { copy->SetBinContent( i+1, h->GetBinContent(i+1)/norm ); }

    auto cv = new TCanvas( "cv1", "cv1", 900, 800 );
    copy->SetMinimum(0.98);

    copy->GetYaxis()->SetTitle("Heartbeat fraction" );
    copy->SetFillStyle(1001);
    copy->SetFillColor(kYellow);
    copy->SetStats(false);

    copy->GetYaxis()->SetTitleOffset(2);
    copy->Draw("hist");

    // add information
    auto text = new TPaveText(0.2, 0.2, 0.5, 0.4, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "Run number: %i", runnumber ) );
    text->Draw();

    gPad->SetLeftMargin(0.2);

    pdfDocument.Add(cv);

  }
}
