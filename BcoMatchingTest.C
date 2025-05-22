#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <micromegas/MicromegasRawDataTimingEvaluation.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void BcoMatchingTest()
{

//   int runnumber = 43817;
//   int nsegments = 42;

//   int runnumber = 43402;
//   int nsegments = 5;

  int runnumber = 45288;
  int nsegments = 10;

  PdfDocument pdfDocument( Form( "Figures/BcoMatchingTest-%08i-0000.pdf", runnumber ));

  gStyle->SetOptStat(0);

  for( int i = 0; i < nsegments; ++i )
  {
    const auto inputfile = Form( "Rootfiles/ScanBCO-%08i-%04i.root", runnumber, i );
    FileManager f( inputfile );
    auto out = f.GetChain( "out" );
    if( !out ) continue;

    const auto hname = Form( "h%i", i );
    auto h = Utils::TreeToHisto( out, hname, "fee_bco_predicted-fee_bco:gtm_bco", "fee_bco>0&&gtm_bco>0", true );
    h->GetXaxis()->SetTitle( "gtm BCO" );
    h->GetYaxis()->SetTitle( "fee BCO (predicted - data)" );
    h->SetTitle(Form( "Run: %05i, segment: %04i", runnumber, i ) );

    auto cvname = Form( "cv_%i", i );
    auto cv = std::make_unique<TCanvas>( cvname, cvname, 800, 800 );
    h->Draw();

    pdfDocument.Add(cv.get());

  }
}
