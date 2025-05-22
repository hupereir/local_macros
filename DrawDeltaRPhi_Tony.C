
const std::string status = "Internal";

//___________________________________________________
void DrawDeltaRPhi_Tony()
{

  set_style( false );

  const TString tag = "_corrected_notpc-new";
  const TString postfix = "";
  const TString inputfile = Form( "Rootfiles/DeltaRPhi_Tony%s%s.root", tag.Data(),postfix.Data());

  auto tf = TFile::Open(inputfile);

  auto cv = new TCanvas( "cv", "cv", 900, 700 );
  auto h = static_cast<TH2*>( tf->Get("h_drphi_2d") );
  h->GetYaxis()->SetTitle( "r#Delta#phi (cm)" );
  h->Draw("col");
  gPad->SetLogz(true);

  auto tg = static_cast<TGraphErrors*>( tf->Get("residuals_layers mean") );
  tg->Draw("P");

  auto text = new TPaveText(0.17,0.73,0.53,0.88, "NDC" );
  text->SetFillColor(0);
  text->SetFillStyle(1001);
  text->SetBorderSize(0);
  text->SetTextAlign(11);
  text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
  text->AddText("p+p #sqrt{s_{NN}} = 200 GeV");
  text->Draw();
}
