
const std::string status = "Internal";

//___________________________________________________
void DrawTrackResiduals()
{

  set_style( false );

  const TString tag = "_corrected";
  const TString postfix = "";
  const TString inputfile = Form( "Rootfiles/TrackResiduals_Tony%s%s.root", tag.Data(),postfix.Data());

  auto tf = TFile::Open( inputfile );

  // x,y profiles
  auto cv = new TCanvas( "cv", "cv", 1400, 700 );
  cv->Divide( 2, 1 );

  const TString labels[] =
  {
    "charge<0",
    "charge>0"
  };

  for( int icut = 0; icut < 2; ++icut )
  {
    auto h0 = static_cast<TH1*>( tf->Get(  Form("hxy_cluster_%i", icut )));
    auto h1 = static_cast<TH1*>( tf->Get(  Form("hxy_state_%i", icut )));

    cv->cd(icut+1);
    h0->SetStats(false);
    h0->Draw();
    h1->Draw("same");

    auto text = new TPaveText(0.17,0.73,0.53,0.88, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
    text->AddText("p+p #sqrt{s_{NN}} = 200 GeV");
    text->AddText( labels[icut] );
    text->Draw();

//     {
//       auto el = new TEllipse(0.,0.,5,5);
//       el->SetFillStyle(0);
//       el->Draw();
//     }
//
//     {
//       auto el = new TEllipse(0.,0.,12,12);
//       el->SetFillStyle(0);
//       el->Draw();
//     }
//
//     {
//       auto el = new TEllipse(0.,0.,75,75);
//       el->SetFillStyle(0);
//       el->Draw();
//     }

  }
}
