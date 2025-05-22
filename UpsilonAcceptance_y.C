#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


double get_rapidity( double e, double pz ) { return 0.5*std::log( (e+pz)/(e-pz) ); }

//____________________________________________________________________________
void UpsilonAcceptance_y()
{
  set_style( false );

  // input
  const TString tag = "_upsilon_pythia8_gen";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*_?.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/UpsilonAcceptance_y%s.pdf", tag.Data() );

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  const TString var( "get_rapidity( _mother._e, _mother._pz )" );
  const TCut base_cut( "_mother._pid == 553&& _daughter[0]._pid == 11 && _daughter[1]._pid == -11");
  const TCut rapidity_cut( "fabs(get_rapidity( _mother._e, _mother._pz ))<=0.5" );

  const TCut trigger_cut[] =
  {
    "",
    rapidity_cut
  };

  const int color[] = { kGray+2,  kBlue };
  const int ncuts = 2;

  std::vector<int> entries;

  const TString label[] =
  {
    "all #Upsilon",
    "|y_{#Upsilon}|<0.5"
  };

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  auto legend = new TLegend( 0.3, 0.75, 0.97, 0.9, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  for( int icut = 0; icut < ncuts; ++icut )
  {
    const auto hname = Form( "h_pt_%i", icut );
    auto h = new TH1F( hname, "", 100, -4, 4 );
    Utils::TreeToHisto( tree, h->GetName(), var, base_cut && trigger_cut[icut], false );

    h->SetTitle("");
    h->GetXaxis()->SetTitle("y");
    h->GetYaxis()->SetTitle("A.U.");
    h->SetMarkerStyle( 20 );
    h->SetMarkerColor(color[icut]);
    h->SetLineColor(color[icut]);
    h->SetLineWidth(2);

    h->SetFillStyle(3003);
    h->SetFillColor(color[icut]);

    if( icut == 0 )
    {
      h->Draw();
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
