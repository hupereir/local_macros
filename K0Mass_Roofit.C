#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooStats/SPlot.h>


R__LOAD_LIBRARY(libRootUtilBase.so)

//______________________________________________
void K0Mass_Roofit()
{

  set_style(kFALSE);
  TGaxis::SetMaxDigits(3);

  bool drawPulls = true;

  const std::vector<int> runnumbers = {53877};
  // const std::vector<int> runnumbers = {53534,53631,53687,53716,53738,53739,53741,53742,53743,53744,53756,53783,53876,53877,53879};
  const auto runLabel = Utils::GetRunLabel(runnumbers);
  const auto runPostfix = Utils::GetRunPostfix(runnumbers);

  const TString tag = "_test1-new";

  const auto pdffilename = drawPulls ?
    Form( "Figures/K0Mass_Roofit%s_Pulls%s.pdf", tag.Data(), runPostfix.Data() ):
    Form( "Figures/K0Mass_Roofit%s%s.pdf", tag.Data(), runPostfix.Data() );

  PdfDocument pdfDocument(pdffilename);

  // load files
  FileManager filemanager;
  for( const auto& runnumber:runnumbers )
  {
    const auto inputFile = Form( "Rootfiles/KFP/k0s%s/KFP-%08i-merged.root", tag.Data(), runnumber );
    filemanager.AddFiles(inputFile);
  }

  // input tree
  auto tree = filemanager.GetChain( "DecayTree" );

  // cuts
  const TCut mvtx_states_cut("min(track_1_MVTX_nStates,track_2_MVTX_nStates)>=2");
  const TCut tpc_states_cut("min(track_1_TPC_nStates,track_2_TPC_nStates)>25");
  const TCut track_chi2_cut(
    "(track_1_nDoF>0&&track_1_chi2/track_1_nDoF>=0&&track_1_chi2/track_1_nDoF<100)&&"
    "(track_2_nDoF>0&&track_2_chi2/track_2_nDoF>=0&&track_2_chi2/track_2_nDoF<100)"
    );
  const TCut track_pt_cut( "min(track_1_pT,track_2_pT)>=0.2");
  const TCut track_ip_cut( "min(fabs(track_1_IP_xy), fabs(track_2_IP_xy))>0.05");
  const TCut track_dca_cut("track_1_track_2_DCA<0.15");
  const TCut decay_length_cut("fabs(K_S0_decayLength)>0.2");

  const TCut cut= mvtx_states_cut && tpc_states_cut&&track_chi2_cut
    &&track_pt_cut
    &&track_ip_cut
    &&track_dca_cut
    &&decay_length_cut
    &&TCut();

  TTree* copy = tree->CopyTree(cut);
  std::cout << "K0Mass_Roofit - done with copy" << std::endl;

  TTree* clone = copy->CloneTree(-1);
  std::cout << "K0Mass_Roofit - done with clone" << std::endl;

  // ROO variable
  const TString var = "K_S0_mass";
  const double mass_min = 0.4;
  const double mass_max = 0.6;
  RooRealVar mass(var, "mass", mass_min, mass_max);

  // mass cut
  TCut mass_cut = Form( "%s>=%f && %s<=%f", var.Data(), mass_min, var.Data(), mass_max );
  std::cout << "K0Mass_Roofit - cut: " << mass_cut << std::endl;

  // data set
  RooDataSet dataSet(var, "data", mass, RooFit::Import(*clone));

  // signal model
  RooRealVar  mean("mean", "mean", 0.498, 0.485, 0.5);
  RooRealVar  sigma("sigma", "sigma", 0.001, 0.000, 0.02);
  RooRealVar  width("width", "Lorentzian Width", 0.003, 0.001, 0.01);

  // left tail
  RooRealVar  cb_alpha_l("alphal", "alphal", 1, 0.000, 5);
  RooRealVar  cb_n_l("nl", "nl", 1.01, 1.0, 20);

  // right tail
  RooRealVar  cb_alpha_r("alphar", "alphar", 1, 0.000, 5);
  RooRealVar  cb_n_r("nr", "nr", 1.01, 1.0, 20);

  // RooGaussian signal("signal", "signal", mass, mean, sigma);
  // RooVoigtian signal("signal", "signal", mass, mean, sigma, width);
  RooCrystalBall signal("signal", "signal", mass, mean, sigma, cb_alpha_l, cb_n_l, cb_alpha_r, cb_n_r);

  // background model
  RooRealVar expConstOne("expConstOne", "expConstOne", -10, -1e2, 0.);
  RooExponential background("background", "background", mass, expConstOne);

  // magnitudes
  RooRealVar  f_signal("f_signal", "f_signal", 1, 0, 1);

  // do the fit
  RooArgList fitModelList(signal, background);
  RooArgList fitFracList(f_signal);
  RooAddPdf model("model", "model", fitModelList, fitFracList);
  model.fitTo(dataSet);

  // get fit result and errors
  const double yield_signal = dataSet.numEntries() * f_signal.getValV();
  const double yield_signal_error = dataSet.numEntries() * f_signal.getError();

  const double mean_value = mean.getValV();
  const double mean_error = mean.getError();

  const double sigma_value = sigma.getValV();
  const double sigma_error = sigma.getError();

  const double yield_background = dataSet.numEntries() * (1.-f_signal.getValV());
  const double yield_background_error = dataSet.numEntries() * f_signal.getError();

  // create untitled RooFit frame for plotting
  auto frame = mass.frame(RooFit::Title(""));

  // data
  RooBinning bins(mass_min, mass_max);
  const int nbins = 50;
  bins.addUniform(nbins, mass_min, mass_max);

  const double bin_width=(mass_max - mass_min)/nbins;

  dataSet.plotOn(frame, RooFit::Binning(bins), RooFit::XErrorSize(0), RooFit::DataError(RooAbsData::SumW2));

  // signal+background
  model.plotOn(
    frame,
    RooFit::Components(RooArgSet(signal, background)),
    RooFit::LineColor(kAzure+8), RooFit::DrawOption("F"),
    RooFit::FillColor(kAzure+8));

  // background
  model.plotOn(frame,
    RooFit::Components(background),
    RooFit::LineColor(kGray),
    RooFit::DrawOption("F"),
    RooFit::FillColor(kGray));

  // signal + background, full line
  model.plotOn(frame, RooFit::LineColor(kBlack));

  // data (again)
  dataSet.plotOn(frame, RooFit::DrawOption("PE1"), RooFit::Binning(bins),RooFit::XErrorSize(0), RooFit::DataError(RooAbsData::SumW2));

  // pulls frame
  auto pulls = frame->pullHist();
  auto frame_pulls = mass.frame(RooFit::Title(""));
  frame_pulls->addPlotable(pulls,"PE1");

  // canvas
  auto cv = new TCanvas( "cv", "", 900, 900 );

  // main pad
  TPad* mainPad = drawPulls ? new TPad("mainPad", "mainPad", 0., 0.3, 1., 1.):cv;
  mainPad->SetTopMargin(0.08);
  mainPad->SetLeftMargin(0.12);

  if( drawPulls )
  {
    mainPad->SetBottomMargin(0);
    mainPad->Draw();
  }

  {
    // main plot
    mainPad->cd();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle( "m(#pi^{+}#pi^{-}) [GeV]");
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitle( Form("candidates/(%.0f MeV)", bin_width*1000) );
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->SetMinimum(0.001);
    frame->Draw();
  }

  // pulls
  if( drawPulls )
  {
    cv->cd();
    auto pullPad = new TPad("pullPad", "pullPad", 0., 0.0, 1., 0.3);
    pullPad->SetTopMargin(0);
    pullPad->SetBottomMargin(0.3);
    pullPad->SetLeftMargin(0.12);
    pullPad->Draw();

    {
      pullPad->cd();
      frame_pulls->SetTitle("");
      frame_pulls->GetXaxis()->SetTitle( "m(#pi^{+}#pi^{-}) [GeV]");
      frame_pulls->GetXaxis()->SetNdivisions(505);
      frame_pulls->GetXaxis()->SetTitleSize(0.12);
      frame_pulls->GetXaxis()->SetLabelSize(0.12);

      frame_pulls->GetYaxis()->SetTitle( "pull" );
      frame_pulls->GetYaxis()->SetNdivisions(505);
      frame_pulls->GetYaxis()->SetTitleSize(0.12);
      frame_pulls->GetYaxis()->SetLabelSize(0.12);
      frame_pulls->GetYaxis()->SetTitleOffset(0.5);
      frame_pulls->GetYaxis()->SetRangeUser(-5.5, 5.5);
      frame_pulls->Draw();
      pullPad->Update();

      Draw::HorizontalLine(pullPad,0)->Draw();
      Draw::HorizontalLine(pullPad,-3)->Draw();
      Draw::HorizontalLine(pullPad,3)->Draw();
    }
  }

  {
    // legend
    mainPad->cd();
    auto legend = new TLegend(0.6,0.66,0.89,0.9);
    legend->AddEntry(frame->findObject("h_K_S0_mass"),"Data","PE2");
    legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]"),"Fit","L");
    legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[signal,background]"),"K^{0}_{S}#rightarrow#pi^{+}#pi^{-}","f");
    legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[background]"),"Comb. Bkg.","f");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->Draw();
  }

  {
    // information
    mainPad->cd();
    auto pave = new TPaveText(0.58,0.29,0.96,0.59, "NDC");
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetTextAlign(11);
//    pave->AddText("#it{#bf{sPHENIX}} Internal");
//    pave->AddText("#it{p}+#it{p} #sqrt{s} = 200 GeV");
    pave->AddText("Early Calibration");
    pave->AddText(Form( "entries = %lli", copy->GetEntries()));
    pave->AddText(Form( "yield = %.0f #pm %.0f", yield_signal, yield_signal_error));
    pave->AddText(Form( "#mu = %.1f #pm %.1f MeV", 1e3*mean_value, 1e3*mean_error));
    pave->AddText(Form( "#sigma = %.1f #pm %.1f MeV", 1e3*sigma_value, 1e3*sigma_error));
    pave->Draw();
  }

  if( false )
  {
    // date
    mainPad->cd();
    auto pave_date = new TPaveText(0.74,0.92,1.,1.,"NDC");
    pave_date->SetFillColor(0);
    pave_date->SetFillStyle(0);
    pave_date->SetBorderSize(0);
    pave_date->SetTextFont(42);
    pave_date->AddText( Utils::GetDate() );
    pave_date->Draw();
  }

  pdfDocument.Add(cv);
}
