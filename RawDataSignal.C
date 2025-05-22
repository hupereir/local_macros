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
TString RawDataSignal(int runNumber = 14091)
{
  
  MicromegasMapping mapping;
  
  const TString inputFile = Form( "DST/MicromegasRawDataEvaluation-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataSignal-%08i-0000.pdf", runNumber );
  const TString rootfileName = Form( "Rootfiles/RawDataSignal-%08i-0000.root", runNumber );
    
  std::cout << "RawDataSignal - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataSignal - pdfFile: " << pdfFile << std::endl;
  std::cout << "RawDataSignal - rootfileName: " << rootfileName << std::endl;

  RootFile rootFile( rootfileName );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const TString var( "samples.adc:samples.strip" );
  const TString var3d = Form( "%s:(samples.tile+8*(samples.layer-55))", var.Data() );
  auto h3d = new TH3I( "h3d", "h3d", 16, 0, 16, 256, 0, 256, 1024, 0, 1024 );
  h3d->GetXaxis()->SetTitle( "tile_id" );
  h3d->GetYaxis()->SetTitle( "strip" );
  h3d->GetZaxis()->SetTitle( "adc" );
  
  Utils::TreeToHisto( tree, h3d->GetName(), var3d, TCut(), false );
  rootFile.Add(h3d);
  
  for( int ilayer = 0; ilayer <2; ++ilayer )
  {
    for( int tile = 0; tile <8; ++tile )
    {
      int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z; 
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, tile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      const auto hname = Form( "h_%i_%i", ilayer, tile );
      const auto htitle = Form( "%s - %i,%i", name.c_str(), ilayer, tile );

      int bin = tile+8*ilayer;
      h3d->GetXaxis()->SetRange( bin+1, bin+1 );
      auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
      h2d->SetName( hname );
      h2d->SetTitle( htitle );
      
      const auto cvname = Form( "cv_%i_%i", ilayer, tile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );
      h2d->Draw("colz");
      pdfDocument.Add( cv );
      
      rootFile.Add(h2d);
    }
  }
  return pdfFile; 
}
