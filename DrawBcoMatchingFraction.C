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


void DrawBcoMatchingFraction( int runnumber=53285)
{
  const TString filename = "/sphenix/data/data02/sphnxpro/production/run3auau/cosmics/new_nocdbtag_v000/DST_STREAMING_EVENT_TPOT/run_00057600_00057700/hist/HIST_DST_STREAMING_EVENT_TPOT_run3auau_new_nocdbtag_v000-00057673-00000.root";
  //  const TString filename = "Fun4All_ReadCombinedData_hp_QA.root";
  std::cout << "DrawBcoMatchingFraction - filename: " << filename << std::endl;

    // create plots
  PdfDocument pdfDocument(Form("Figures/BcoMacthingFraction-%08i-0000.pdf", runnumber));

  auto tfile = TFile::Open(filename);

  if( true )
  {
    auto h = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_npacket_bco" ) );
    auto cv = new TCanvas( "cv0", "cv0", 900, 800 );
    h->Draw("hist");
    gPad->SetLogy(true);

    // calculate mean
    double sum = 0;
    double count = 0;
    for( int i = 0 ; i < h->GetNbinsX(); ++i )
    {
      sum += double(i) * h->GetBinContent(i+1);
      count +=  h->GetBinContent( i+1);
    }

    const double mean = sum/count;

    // add information
    auto text = new TPaveText(0.6, 0.6, 0.9, 0.8, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "Run number: %i", runnumber ) );
    text->Draw();

    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto h = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_packet_stat" ) );
    auto copy = get_copy(h);

    // normalize
    auto norm = h->GetBinContent(1);

    std::cout << "DrawPackets - triggers: " << norm << std::endl;

    for( int i = 0; i < 4; ++i )
    { copy->SetBinContent( i+1, h->GetBinContent(i+1)/norm ); }

    auto cv = new TCanvas( "cv1", "cv1", 900, 800 );
    copy->SetMinimum(0.);

    copy->GetYaxis()->SetTitle("GL1 trigger fraction" );
    copy->SetFillStyle(1001);
    copy->SetFillColor(kYellow);
    copy->SetStats(false);

    copy->GetYaxis()->SetTitleOffset(1.4);
    copy->Draw("hist");

    // add information
    auto text = new TPaveText(0.2, 0.2, 0.5, 0.4, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "Run number: %i", runnumber ) );
    text->AddText( Form( "Triggers: %.3g", double(norm) ) );
    text->AddText( Form( "fraction 5001: %.4g %%", 100.*double(h->GetBinContent(2)) ) );
    text->AddText( Form( "fraction 5002: %.4g %%", 100.*double(h->GetBinContent(3)) ) );
    text->AddText( Form( "fraction all: %.4g %%", 100.*double(h->GetBinContent(4)) ) );
    text->Draw();

    gPad->SetLeftMargin(0.16);

    pdfDocument.Add(cv);

  }

  if( true )
  {
    auto htot = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_waveform_count_total" ) );
    auto hdropped = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_waveform_count_dropped_bco" ) );
    auto copy = get_copy(hdropped);

    for( int i = 0; i < htot->GetNbinsX(); ++i )
    {
      copy->SetBinContent( i+1, double(hdropped->GetBinContent(i+1))/htot->GetBinContent(i+1) );
    }

    copy->SetTitle( "dropped waveform fraction per packet" );
    copy->GetXaxis()->SetTitle( "packet id" );
    copy->GetYaxis()->SetTitle( "dropped waveform fraction (bco mismatch)" );

    auto cv = new TCanvas( "cv2", "cv2", 900, 800 );
    copy->SetMinimum(0.);
    copy->SetMaximum(1.0);

    copy->SetFillStyle(1001);
    copy->SetFillColor(kYellow);
    copy->SetStats(false);

    copy->GetYaxis()->SetTitleOffset(1.4);
    copy->Draw("hist");

    // add information
    auto text = new TPaveText(0.2, 0.2, 0.5, 0.4, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "Run number: %i", runnumber ) );
    text->AddText( Form( "fraction 5001: %.4g %%", 100.*double(copy->GetBinContent(1)) ) );
    text->AddText( Form( "fraction 5002: %.4g %%", 100.*double(copy->GetBinContent(2)) ) );
    text->Draw();

    gPad->SetLeftMargin(0.16);

    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto htot = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_waveform_count_total" ) );
    auto hdropped = static_cast<TH1*>( tfile->Get("h_MicromegasBCOQA_waveform_count_dropped_pool" ) );
    auto copy = get_copy(hdropped);

    for( int i = 0; i < htot->GetNbinsX(); ++i )
    {
      copy->SetBinContent( i+1, double(hdropped->GetBinContent(i+1))/htot->GetBinContent(i+1) );
    }

    copy->SetTitle( "dropped waveform fraction per packet" );
    copy->GetXaxis()->SetTitle( "packet id" );
    copy->GetYaxis()->SetTitle( "dropped waveform fraction (pool mismatch)" );

    auto cv = new TCanvas( "cv3", "cv3", 900, 800 );
    copy->SetMinimum(0.);
    copy->SetMaximum(0.5);

    copy->SetFillStyle(1001);
    copy->SetFillColor(kYellow);
    copy->SetStats(false);

    copy->GetYaxis()->SetTitleOffset(1.4);
    copy->Draw("hist");

    // add information
    auto text = new TPaveText(0.2, 0.2, 0.5, 0.4, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "Run number: %i", runnumber ) );
    text->AddText( Form( "fraction 5001: %.4g %%", 100.*double(copy->GetBinContent(1)) ) );
    text->AddText( Form( "fraction 5002: %.4g %%", 100.*double(copy->GetBinContent(2)) ) );
    text->Draw();

    gPad->SetLeftMargin(0.16);
    pdfDocument.Add(cv);
  }

}
