#include <RootUtil/PdfDocument.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//_____________________________________________________________________________
TGraphErrors* combine_histograms( TH1* efficiency, TH1* rejection )
{
  auto out = new TGraphErrors;
  int ipt = 0;
  for( int ibin = 0; ibin < efficiency->GetNbinsX(); ++ibin )
  {

    // get values
    const auto eff = efficiency->GetBinContent(ibin+1);
    const auto eff_error = efficiency->GetBinError(ibin+1);
    const auto rej = rejection->GetBinContent(ibin+1);
    const auto rej_error = rejection->GetBinError(ibin+1);

    if(std::isnan(eff_error)) continue;
    if(std::isinf(rej)) continue;

    std::cout << "combine_histograms - (" << eff << "," << eff_error << ") (" << rej << "," << rej_error << ")" << std::endl;

    // assign
    out->SetPoint(ipt, eff, rej);
    out->SetPointError(ipt, eff_error, rej_error);
    ++ipt;
  }

  out->GetXaxis()->SetTitle(efficiency->GetXaxis()->GetTitle());
  out->GetYaxis()->SetTitle(rejection->GetXaxis()->GetTitle());
  return out;
}

//_____________________________________________________________________________
void DrawEOverPCuts()
{

  set_style(false);
  const TString global_tag = "_lpt_p2";

  const TString efficiencyFile = Form("Rootfiles/e_over_p_efficiency%s.root", global_tag.Data());
  const TString rejectionFile = Form("Rootfiles/e_over_p_rejection%s.root", global_tag.Data());

  const TString pdfFile = Form("Figures/e_over_p_cuts%s.pdf", global_tag.Data());
  PdfDocument pdfDocument(pdfFile);

  auto tfile_efficiency = TFile::Open(efficiencyFile);
  auto efficiency_single_electron = static_cast<TH1*>(tfile_efficiency->Get("e_over_p_single_electron_efficiency"));
  auto efficiency_single_positron = static_cast<TH1*>(tfile_efficiency->Get("e_over_p_single_positron_efficiency"));

  auto tfile_rejection = TFile::Open(rejectionFile);
  auto rejection_single_piminus = static_cast<TH1*>(tfile_rejection->Get("e_over_p_single_piminus_rejection"));
  auto rejection_single_piplus = static_cast<TH1*>(tfile_rejection->Get("e_over_p_single_piplus_rejection"));

  static constexpr double max_rejection = 300;

  if( true )
  {
    // electron vs piminus
    auto combined = combine_histograms(efficiency_single_electron,rejection_single_piminus);
    combined->GetXaxis()->SetTitle( "e^{-} efficiency" );
    combined->GetYaxis()->SetTitle( "#pi^{-} rejection" );
    combined->SetMarkerStyle(20);
    combined->SetMarkerColor(1);
    combined->SetLineColor(1);

    combined->GetXaxis()->SetRangeUser(0.5,1);
    combined->SetMaximum(max_rejection);

    // canvas
    auto cv( new TCanvas( "cv", "cv", 900, 900 ) );
    combined->Draw("AP");
    pdfDocument.Add(cv);
  }

  if( true )
  {
    // positron vs piplus
    auto combined = combine_histograms(efficiency_single_positron,rejection_single_piplus);
    combined->GetXaxis()->SetTitle( "e^{+} efficiency" );
    combined->GetYaxis()->SetTitle( "#pi^{+} rejection" );
    combined->SetMarkerStyle(20);
    combined->SetMarkerColor(1);
    combined->SetLineColor(1);

    combined->GetXaxis()->SetRangeUser(0.5,1);
    combined->SetMaximum(max_rejection);

    // canvas
    auto cv( new TCanvas( "cv2", "cv", 900, 900 ) );
    combined->Draw( "AP" );
    pdfDocument.Add(cv);
  }
}
