#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//________________________________________________
TString Occupancy()
{
  const TString inputFile =  "DST/CONDOR_MicromegasEvaluation_hijing/dst_simeval_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000007-000*.root";
  const TString pdfFile = "Figures/Occupancy_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000007.pdf";
  const TString rootfilename = "Rootfiles/Occupancy_sHijing_0_20fm_bkg_0_20fm-0000000007.root";

//   const TString inputFile =  "DST/CONDOR_MicromegasEvaluation_hijing/dst_simeval_sHijing_0_20fm-0000000007-000*.root";
//   const TString pdfFile = "Figures/Occupancy_sHijing_0_20fm-0000000007.pdf";
//   const TString rootfilename = "Rootfiles/Occupancy_sHijing_0_20fm-0000000007.root";

  PdfDocument pdfDocument( pdfFile );
  RootFile rootfile( rootfilename );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const TString var_raw( "_tiles[]._strips_total" );
  const TString var(  "100.*_tiles[]._strips_total/256." );
  const TString var_2d = Form( "%s:100*pow(_events[0]._bimp/12.8,2)", var.Data() );
  const TCut cut[] = {
    "_tiles[]._layer == 55",
    "_tiles[]._layer == 56"
  };

  const TString label[] = {
    "occupancy (#Phi strips) (%)",
    "occupancy (#it{Z} strips) (%)"
  };

  static const int ncut =  2;
  for( int i =0; i<ncut; ++i )
  {

    {
      auto hname = Form( "h_%i", i );
      auto h = new TH1F( hname, "", 256, 0, 256 );
      h->GetXaxis()->SetTitle( "hit count" );
      Utils::TreeToHisto( tree, h->GetName(), var_raw, cut[i], false );
      rootfile.Add( h );
      std::cout << "Occupancy - mean: " << 100.*h->GetMean()/256 << std::endl;
    }

    auto hname = Form( "h2_%i", i );
    auto h = new TH2F( hname, "", 10, 0, 100, 100, 0, 100 );
    h->GetXaxis()->SetTitle( "centrality (%)" );
    h->GetYaxis()->SetTitle( label[i] );
    Utils::TreeToHisto( tree, h->GetName(), var_2d, cut[i], false );
    rootfile.Add( h );

    auto pname = Form( "p_%i", i );
    auto p = new TProfile( pname, "", 10, 0, 100 );
    p->GetXaxis()->SetTitle( "centrality (%)" );
    p->GetYaxis()->SetTitle( label[i] );
    p->SetMarkerStyle( 20 );
    Utils::TreeToHisto( tree, p->GetName(), var_2d, cut[i], false );
    rootfile.Add( p );

    // create canvas and divide
    auto cv = new TCanvas( "cv", "cv", 980, 900 );
    gPad->SetTopMargin( 0.05 );
    gPad->SetRightMargin( 0.05);

    h->Draw( "col" );
    p->Draw( "same" );
    pdfDocument.Add(cv);
  }
  return pdfFile;
}
