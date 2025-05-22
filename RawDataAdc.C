#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

//_____________________________________________________________________________
TString RawDataAdc( int runNumber = 7372 )
{

  MicromegasMapping mapping;

  const TString inputFile = Form( "DST/MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataAdc-%08i-0000.pdf", runNumber );
  const TString rootfileName = Form( "Rootfiles/RawDataAdc-%08i-0000.root", runNumber );

  std::cout << "RawDataTiming - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataTiming - pdfFile: " << pdfFile << std::endl;
  std::cout << "RawDataTiming - rootfileName: " << rootfileName << std::endl;

  RootFile rootFile( rootfileName );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const TString var( "adc:sample" );
  const TString var3d = Form( "%s:((strip/64)+4*(tile+8*(layer-55)))", var.Data() );
  auto h3d = new TH3I( "h3d", "h3d", 64, 0, 64, 300, 0, 300, 1100, 0, 1100 );
  h3d->GetXaxis()->SetTitle( "region_id" );
  h3d->GetYaxis()->SetTitle( "sample" );
  h3d->GetZaxis()->SetTitle( "adc" );
  Utils::TreeToHisto( tree, h3d->GetName(), var3d, TCut(), false );
  rootFile.Add(h3d);


  for( int ilayer = 0; ilayer <2; ++ilayer )
  {
    for( int itile = 0; itile <8; ++itile )
    {
      const auto cvname = Form( "cv_%i_%i", ilayer, itile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );
      cv->Divide(2,2);

      for( int iregion = 0; iregion < 4; ++iregion )
      {

        const int layer = 55+ilayer;
        const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
        const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
        const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );

        // region are named opposite to strips.
        const int region = 4-iregion;
        const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
        const auto htitle = Form( "%s_R%i - %i,%i", name.c_str(), region, layer, itile );

        int bin = iregion + 4*(itile + 8*ilayer );
        h3d->GetXaxis()->SetRange( bin+1, bin+1 );
        auto h1d = static_cast<TH1*>( h3d->Project3D( "z" ) );
        h1d->SetName( hname );
        h1d->SetTitle( htitle );

        cv->cd(iregion+1);
        h1d->Draw();
        gPad->SetLogy(true);
        rootFile.Add(h1d);
      }

      pdfDocument.Add( cv );

    }
  }

  return pdfFile;

}
