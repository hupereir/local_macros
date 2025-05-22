#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
#include <RootUtil/Color.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
void DrawJPsiAcceptance_new()
{
  set_style( false );

  // input
  const TString tag = "_jpsi_pythia8_gen";
  const TString inputFile = Form( "Rootfiles/JpsiAcceptance_new%s.root", tag.Data() );

  auto tfile = TFile::Open( inputFile );

  const int color[] = { kGray+2,  kBlue, kGreen+2, kOrange+2, kRed };
  const int ncuts = 5;

  std::vector<int> entries;

  const TString label[] =
  {
    "all J/#psi",
    "J/#psi in acceptance",
    "trigger selection",
    "electron ID (E>1 GeV)",
    "electron ID (E>2 GeV)"
  };

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  cv->SetLeftMargin( 0.13 );

  auto legend = new TLegend( 0.46, 0.72, 0.96, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  bool first = true;
  for( int icut = 1; icut < ncuts; ++icut )
  {
    const auto hname = Form( "h_pt_%i", icut );
    auto h = static_cast<TH1*>( tfile->Get(hname) );

    h->GetXaxis()->SetTitleOffset(1.3) ;

    h->SetTitle("");
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("A.U.");
    h->SetMarkerStyle( 20 );
    h->SetMarkerColor(color[icut]);

    h->SetLineColor(color[icut]);
    h->SetLineWidth(2);

    h->SetFillStyle(1001);
    const int merge = Color( color[icut] ).Merge( kWhite, 0.1 );
    h->SetFillColor(merge);

    if( first )
    {
      first = false;
      h->Draw();
      gPad->SetLogy( true );
      legend->Draw();
    } else h->Draw("same");

    legend->AddEntry( h, label[icut], "L" );
    entries.push_back( h->GetEntries() );

  }

  auto arrow = new TArrow(0.3,0.625, 0.3, 0.7, 0.03, "|>" );
  arrow->SetNDC( true );
  arrow->SetLineWidth(4);
  arrow->Draw();

  {
    auto text = new TPaveText(0.35,0.63,0.53,0.69, "NDC" );
    text->SetFillColorAlpha(0,0.8);
    text->SetFillStyle(1001);
    text->SetBorderSize(0);
    text->AddText("#times 5 (ML)");
    text->Draw();
  }

  cv->SaveAs("Figures/JpsiAcceptance_new.pdf" );

  for( size_t i = 0; i < entries.size(); ++i )
  {
    std::cout << "JpsiAcceptance - i=" << i << " entries: " << entries[i] << " ratio: " << double(entries[i])/entries[0] << std::endl;
  }

}
