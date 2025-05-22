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
R__LOAD_LIBRARY(libg4eval_hp.so)

//______________________________________________
void JPsiMass_Roofit()
{

  set_style(kFALSE);
  TGaxis::SetMaxDigits(3);

  bool drawPulls = true;

  // load data
  const TString tag = "_jpsi_pythia8";

  const auto pdffilename = drawPulls ?
    Form( "Figures/JPsiMass_Roofit%s_Pulls.pdf", tag.Data()):
    Form( "Figures/JPsiMass_Roofit%s.pdf", tag.Data());

  PdfDocument pdfDocument(pdffilename);

  // file manager
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*_??.root", tag.Data(), tag.Data() );
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return;

  const TCut pattern_cut(
    "min(_track_pairs._tracks[0]._nclusters_mvtx,_track_pairs._tracks[1]._nclusters_mvtx)>2"
    "&&min(_track_pairs._tracks[0]._nclusters_tpc,_track_pairs._tracks[1]._nclusters_tpc)>30"
    );

  const TCut track_pt_cut( "min(_track_pairs._tracks[0]._pt, _track_pairs._tracks[1]._pt)>=0.2" );
  const TCut cut = pattern_cut&&track_pt_cut;

  auto copy = tree->CopyTree(cut);
  std::cout << "JPsiMass_Roofit - done with copy - entries: " << copy->GetEntries() << std::endl;

  TrackingEvaluator_hp::Container* container = nullptr;
  copy->SetBranchAddress("DST#EVAL#TrackingEvaluator_hp::Container", &container );

  // create a simple tree with only the mass, and project manually
  auto clone = new TTree("T2", "T2");
  {
    double mass = 0;
    clone->Branch("mass", &mass );
    for( int i=0; i<copy->GetEntries(); ++i )
    {
      copy->GetEntry(i);
      for( const auto& trackpair:container->trackPairs() )
      {
        mass = trackpair._m;
        clone->Fill();
      }
    }
  }

  // ROO variable
  const TString var = "mass";
  const double mass_min = 2.3;
  const double mass_max = 3.6;
  RooRealVar mass(var, "mass", mass_min, mass_max);

  // mass cut
  TCut mass_cut = Form( "%s>=%f && %s<=%f", var.Data(), mass_min, var.Data(), mass_max );
  std::cout << "JPsiMass_Roofit - cut: " << mass_cut << std::endl;

  // data set
  RooDataSet dataSet(var, "data", mass, RooFit::Import(*clone));

  // signal model
  static constexpr float mass_jpsi = 3.09 ;
  RooRealVar  mean("mean", "mean", mass_jpsi, mass_jpsi-0.2, mass_jpsi+0.2);
  RooRealVar  sigma("sigma", "sigma", 0.05, 0.000, 0.1);
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
    frame->GetXaxis()->SetTitle( "m(e^{+}e^{-}) [GeV]");
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
      frame_pulls->GetXaxis()->SetTitle( "m(e^{+}e^{-}) [GeV]");
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
    auto legend = new TLegend(0.68,0.63,0.97,0.87);
    legend->AddEntry(frame->findObject("h_mass"),"Data","PE2");
    legend->AddEntry(frame->findObject("model_Norm[mass]"),"Fit","L");
    legend->AddEntry(frame->findObject("model_Norm[mass]_Comp[signal,background]"),"J/#psi#rightarrowe^{+}e^{-}","f");
    legend->AddEntry(frame->findObject("model_Norm[mass]_Comp[background]"),"Comb. Bkg.","f");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->Draw();
  }

  {
    // information
    mainPad->cd();
    auto pave = new TPaveText(0.13,0.03,0.51,0.33, "NDC");
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetTextAlign(11);
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

