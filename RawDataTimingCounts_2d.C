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
TString RawDataTimingCounts_2d( int runNumber = 20464 )
{

  MicromegasMapping mapping;

  const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/MicromegasCombinedDataEvaluation_masked-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataTimingCounts_2d-%08i-0000.pdf", runNumber );
  const TString rootfileName = Form( "Rootfiles/RawDataTimingCounts_2d-%08i-0000.root", runNumber );

  std::cout << "RawDataTiming - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataTiming - pdfFile: " << pdfFile << std::endl;
  std::cout << "RawDataTiming - rootfileName: " << rootfileName << std::endl;

  RootFile rootFile( rootfileName );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  Utils::max_entries = 2000;

  std::cout << "RawDataTimingCounts_2d - entries: " << tree->GetEntries() << std::endl;

  // project adc vs sample vs strip in a 3D histogram
  // strips are grouped by chuncks of 64 corresponding to the 4 Resist region in each of the detectors
  const TString var( "samples.sample:samples.strip" );
  const TString var3d = Form( "%s:(samples.tile+8*(samples.layer-55))", var.Data() );
  auto h3d = new TH3I( "h3d", "h3d", 16, 0, 16, 256, 0, 256, 150, 0, 150 );
  h3d->GetXaxis()->SetTitle( "region_id" );
  h3d->GetYaxis()->SetTitle( "strip" );
  h3d->GetZaxis()->SetTitle( "sample" );

  TCut cut( "samples.adc>samples.pedestal+5.*samples.rms" );
  Utils::TreeToHisto( tree, h3d->GetName(), var3d, cut, false );
  rootFile.Add(h3d);

  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer <2; ++ilayer )
  {

    // get actual layer and segmentation
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;

    for( int itile = 0; itile <8; ++itile )
    {
      // generate hitset key from layer, tile and segmentation. This allows to retrieve the detector name from the Mapping
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );

      // create canvas and divide
      const auto cvname = Form( "cv_%i_%i", ilayer, itile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );

      // due to mapping internals, region are named opposite to strips.
      const auto hname = Form( "h_%i_%i", ilayer, itile );
      const auto htitle = Form( "%s - %i,%i", name.c_str(), layer, itile );

      // get the bin matching layer, tile and region in the 3D histgoram
      const int bin = itile + 8*ilayer;
      h3d->GetXaxis()->SetRange( bin+1, bin+1 );

      // get the corresponding 2D histogram (adc vs sample)
      auto h = static_cast<TH1*>( h3d->Project3D("zy") );
      h->SetName( hname );
      h->SetTitle( htitle );

      // draw
      h->Draw( "colz" );

      gPad->SetRightMargin( 0.2 );

      rootFile.Add(h);

      // add to pdfdocument
      pdfDocument.Add( cv );

    }
  }

  return pdfFile;

}


//_______________________________________________________
void process_all()
{
  using run_map_t = std::map<int,int>;
  run_map_t run_map =
  {
    { 50, 20464 },
    { 100, 20465 },
    { 150, 20466 },
    { 200, 20467 },
    { 250, 20468 },
    { 300, 20469 },
    { 350, 20470 },
    { 400, 20471 }
  };

  for( const auto& [hv, runnumber]:run_map )
  { RawDataTimingCounts_2d(runnumber); }
}
