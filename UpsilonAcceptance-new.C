#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
void UpsilonAcceptance_new()
{
  set_style( false );

  // input
  const TString tag = "_upsilon_pythia8_gen";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*_?.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/UpsilonAcceptance_new%s.pdf", tag.Data() );

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  const TString var( "_mother._pt" );
  const TCut base_cut( "_mother._pid == 553&& _daughter[0]._pid == 11 && _daughter[1]._pid == -11");
  const TCut acceptance_cut( "fabs(_daughter[0]._eta)<1 && fabs(_daughter[1]._eta)<1" );
  const TCut trigger_cut( "_daughter[0]._e>4 || _daughter[1]._e>4" );
  const TCut eid_cut( "_daughter[0]._e>2 && _daughter[1]._e>2" );
  // const TCut eid_cut( "_daughter[0]._e>1 && _daughter[1]._e>1" );

  const TCut selection_cut[] =
  {
    "",
    acceptance_cut,
    acceptance_cut && trigger_cut,
    acceptance_cut && trigger_cut&&eid_cut
  };

  const int color[] = { kGray+2,  kBlue, kGreen+2, kRed };
  const int ncuts = 4;

  std::vector<int> entries;

  const TString label[] =
  {
    "all #Upsilon",
    "#Upsilon in acceptance (both electron |#eta|>1)",
    "trigger selection (one electron E>4 GeV)",
    "electron identification (both electrons E>2 GeV)"
    // "electron identification (both electrons E>1 GeV)"
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
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("A.U.");
    h->SetMarkerStyle( 20 );
    h->SetMarkerColor(color[icut]);
    h->SetLineColor(color[icut]);
    h->SetLineWidth(2);

    h->SetFillStyle(3003);
    h->SetFillColor(color[icut]);

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

  pdfDocument.Add( cv );

  for( size_t i = 0; i < entries.size(); ++i )
  {
    std::cout << "UpsilonAcceptance - i=" << i << " entries: " << entries[i] << " ratio: " << double(entries[i])/entries[0] << std::endl;
  }

}
