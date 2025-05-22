#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
void GtmBCO( int runnumber = 45490, int segment = 0 )
{
  gStyle->SetOptStat(0);

  const TString inputfile = Form( "Rootfiles/ScanBCO-%08i-%04i.root", runnumber, segment );
  const TString pdfFile = Form( "Figures/GtmBCO-%08i-%04i.pdf", runnumber, segment );

//   const TString inputfile( "Rootfiles/ScanBCO-00043817-0000.root" );
//   const TString pdfFile( "Figures/GtmBCO-00043817-0000.pdf" );

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputfile );
  auto tree = fileManager.GetChain( "out" );
  if( !tree ) return TString();

  std::array<int,2> packets = {5001,5002};
  std::array<int,2> color = {1,2};


  bool first= true;
  auto cv = new TCanvas( "cv", "cv", 900, 800 );
  cv->SetTopMargin(0.05);
  cv->SetRightMargin(0.15);

  auto legend = new TLegend( 0.5, 0.84, 0.97, 0.94, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  const TString var = "gtm_bco";
  auto h0 = Utils::TreeToHisto( tree, "h0", var, "gtm_bco>0", true );


  std::array<TH1*,2> harray;

  for( int i = 0; i<2; ++i )
  {
    const auto packet = packets[i];
    const TCut cut = Form( "gtm_bco > 0 && packet_id == %i", packet );

    auto hname = Form( "h_%i", packet );
    auto h = static_cast<TH1*>(h0->Clone( hname ) );
    h->Reset();
    Utils::TreeToHisto( tree, hname, var, cut, false );
    h->SetLineColor( color[i] );
    h->SetTitle("");
    h->GetXaxis()->SetTitle( "GTM BCO" );
    std::cout << "packet: " << packet << " entries: " << h->GetEntries() << std::endl;
    if( first )
    {
      h->Draw();
      first = false;
    } else h->Draw("same");
    harray[i]=h;

    legend->AddEntry( h, Form( "packet: %i", packet ), "L" );
    std::cout << "GtmBCO - packet: " << packet<< " entries: " << h->GetEntries() << std::endl;
  }
  legend->Draw();

  auto max = 1.2*std::max( harray[0]->GetMaximum(), harray[1]->GetMaximum() );
  harray[0]->SetMaximum( max );
  harray[1]->SetMaximum( max );
  cv->Update();
  pdfDocument.Add(cv );
}
