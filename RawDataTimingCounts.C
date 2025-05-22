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

#include <micromegas/MicromegasHotChannelMapData.h>
#include <micromegas/MicromegasMapping.h>

namespace
{
  MicromegasHotChannelMapData hot_channels;
}

//_____________________________________________________________________________
bool is_good( int layer, int tile, int strip )
{ return !hot_channels.is_hot_channel( layer, tile, strip ); }

//_____________________________________________________________________________
TString RawDataTimingCounts( int runNumber = 9467 )
{

    MicromegasMapping mapping;

    static const bool use_mask = true;

    const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/MicromegasCombinedDataEvaluation-%08i-0000-test.root", runNumber );
    // const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/MicromegasCombinedDataEvaluation-%08i-0000-full.root", runNumber );
    const TString pdfFile = use_mask ?
        Form( "Figures/RawDataTimingCounts_masked-%08i-0000.pdf", runNumber ):
        Form( "Figures/RawDataTimingCounts-%08i-0000.pdf", runNumber );

    const TString rootfileName = use_mask ?
        Form( "Rootfiles/RawDataTimingCounts_masked-%08i-0000.root", runNumber ):
        Form( "Rootfiles/RawDataTimingCounts-%08i-0000.root", runNumber );

    std::cout << "RawDataTiming - inputFile: " << inputFile << std::endl;
    std::cout << "RawDataTiming - pdfFile: " << pdfFile << std::endl;
    std::cout << "RawDataTiming - rootfileName: " << rootfileName << std::endl;

    // read hot channels
    if( use_mask )
    {
      const TString calibrationFile = Form( "Calibrations/TPOT_HotChannels-%08i-0000.root", runNumber );
      std::cout << "RawDataTiming - calibrationFile: " << calibrationFile << std::endl;
      hot_channels.read( calibrationFile.Data() );
    }
    RootFile rootFile( rootfileName );
    PdfDocument pdfDocument( pdfFile );

    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "T" );

    // project adc vs sample vs strip in a 3D histogram
    // strips are grouped by chuncks of 64 corresponding to the 4 Resist region in each of the detectors
    const TString var( "samples.sample" );
    const TString var2d = Form( "%s:((samples.strip/64)+4*(samples.tile+8*(samples.layer-55)))", var.Data() );
    auto h2d = new TH2I( "h2d", "h2d", 64, 0, 64, 150, 0, 150 );
    h2d->GetXaxis()->SetTitle( "region_id" );
    h2d->GetYaxis()->SetTitle( "sample" );
    h2d->GetZaxis()->SetTitle( "counts" );

    TCut cut( "samples.adc>samples.pedestal+5.*samples.rms" );

    if( use_mask )
    {
        // Utils::max_entries = 200;
        const TCut hot_channel_cut = "is_good( samples.layer, samples.tile, samples.strip )";
        Utils::TreeToHisto( tree, h2d->GetName(), var2d, cut&&hot_channel_cut, false );
    } else {
        // Utils::max_entries = 200;
        Utils::TreeToHisto( tree, h2d->GetName(), var2d, cut, false );
    }

    rootFile.Add(h2d);

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
            cv->Divide(2,2);

            for( int iregion = 0; iregion < 4; ++iregion )
            {
                // due to mapping internals, region are named opposite to strips.
                const int region = 4-iregion;
                const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
                const auto htitle = Form( "%s_R%i - %i,%i", name.c_str(), region, layer, itile );

                // get the bin matching layer, tile and region in the 3D histgoram
                const int bin = iregion + 4*(itile + 8*ilayer );
                h2d->GetXaxis()->SetRange( bin+1, bin+1 );

                // get the corresponding 2D histogram (adc vs sample)
                auto h = static_cast<TH1*>( h2d->ProjectionY() );
                h->SetName( hname );
                h->SetTitle( htitle );

                // draw
                cv->cd(iregion+1);
                h->Draw();
                rootFile.Add(h);
            }

            // add to pdfdocument
            pdfDocument.Add( cv );

        }
    }

    return pdfFile;

}

//___________________________________
void process_all()
{
  for( const auto& runnumber: { 9438, 20467, 20468, 20469, 20470, 20471, 21082, 21923, 21928, 21933, 21936, 21939, 21950 } )
  { RawDataTimingCounts( runnumber ); }
}
