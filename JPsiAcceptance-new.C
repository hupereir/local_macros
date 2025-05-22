#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>
#include <RootUtil/Color.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString JPsiAcceptance_new()
{
  set_style( false );

  gStyle->SetOptStat(0);

  // input
  const TString tag = "_jpsi_pythia8_gen";
  // const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*_000?.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/JpsiAcceptance_new%s.pdf", tag.Data() );
  const TString rootFilename = Form( "Rootfiles/JpsiAcceptance_new%s.root", tag.Data() );

  PdfDocument pdfDocument( pdfFile );
  RootFile rootfile( rootFilename );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  const TString var( "_mother._pt" );
  const TCut base_cut( "_mother._pid==443&& _daughter[0]._pid == 11 && _daughter[1]._pid == -11");
  const TCut acceptance_cut( "fabs(_daughter[0]._eta)<1 && fabs(_daughter[1]._eta)<1" );
  const TCut trigger_cut( "_daughter[0]._e>4 || _daughter[1]._e>4" );
  const TCut eid_cut_1( "_daughter[0]._e>1 && _daughter[1]._e>1" );
  const TCut eid_cut_2( "_daughter[0]._e>2 && _daughter[1]._e>2" );

  const TCut selection_cut[] =
  {
    "",
    acceptance_cut,
    acceptance_cut && trigger_cut,
    acceptance_cut && trigger_cut&&eid_cut_1,
    acceptance_cut && trigger_cut&&eid_cut_2
  };

  const int color[] = { kGray+2,  kBlue, kGreen+2, kOrange+2, kRed };
  const int ncuts = 5;

  std::vector<int> entries;

  const TString label[] =
  {
    "all J/#psi",
    "J/#psi in acceptance (both #eta<1)",
    "trigger selection (single E>4 GeV)",
    "electron ID (both E>1 GeV)",
    "electron ID (both E>2 GeV)"
  };

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  auto legend = new TLegend( 0.3, 0.75, 0.97, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  bool first = true;
  for( int icut = 1; icut < ncuts; ++icut )
  {
    const auto hname = Form( "h_pt_%i", icut );
    auto h = new TH1F( hname, "", 100, 0, 20 );
    Utils::TreeToHisto( tree, h->GetName(), var, base_cut && selection_cut[icut], false );

    h->SetTitle("");
    h->GetXaxis()->SetTitle("J/#psi #it{p}_{T} (GeV/#it{c})");
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

    rootfile.Add( h );

  }

  if( false )
  {
    auto arrow = new TArrow(0.28,0.52, 0.28, 0.64, 0.03, "|>" );
    arrow->SetNDC( true );
    arrow->SetLineWidth(4);
    // arrow->PaintArrow(4,100, 4, 500, 0.03, "|>" );
    arrow->Draw();

    auto text = new TPaveText(0.3,0.55,0.47,0.61, "NDC" );
    text->SetFillColorAlpha(0,0.8);
    text->SetFillStyle(1001);
    text->SetBorderSize(0);
    text->AddText("#times 5 (ML)");
    text->Draw();
  }

  cv->SaveAs("Figures/JpsiAcceptance_new.pdf" );

  pdfDocument.Add( cv );

  for( size_t i = 0; i < entries.size(); ++i )
  {
    std::cout << "JpsiAcceptance - i=" << i << " entries: " << entries[i] << " ratio: " << double(entries[i])/entries[0] << std::endl;
  }

  return pdfFile;
}
