#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//_____________________________________________________________________________
TString DrawRawDataTimingCounts_2d( int runNumber = 9455 )
{
  const TString rootfileName = Form( "Rootfiles/RawDataTimingCounts_2d-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataTimingCounts_2d-%08i-0000.pdf", runNumber );
  
  auto tfile = TFile::Open( rootfileName, "READ");

  PdfDocument pdfDocument( pdfFile );
  
  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer <2; ++ilayer )
  {
        
    for( int itile = 0; itile <8; ++itile )
    {
      const auto hname = Form( "h_%i_%i", ilayer, itile );
      const auto h = static_cast<TH2*>( tfile->Get( hname ) );

      // create canvas and divide
      const auto cvname = Form( "cv_%i_%i", ilayer, itile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );

      // draw
      h->SetMaximum(500);
      h->Draw( "colz" );
      gPad->SetRightMargin( 0.2 );
      pdfDocument.Add( cv );

    }
  }

  return pdfFile;

}
