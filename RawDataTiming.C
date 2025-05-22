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
TString RawDataTiming( int runNumber = 34726 )
{

    MicromegasMapping mapping;

    static const bool use_mask = true;

    const TString inputFile = Form( "DST/MicromegasRawDataEvaluation-%08i-0000-test.root", runNumber );
    // const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/MicromegasCombinedDataEvaluation-%08i-0000-full.root", runNumber );
    const TString pdfFile = use_mask ?
        Form( "Figures/RawDataTiming_masked-%08i-0000.pdf", runNumber ):
        Form( "Figures/RawDataTiming-%08i-0000.pdf", runNumber );

    const TString rootfileName = use_mask ?
        Form( "Rootfiles/RawDataTiming_masked-%08i-0000.root", runNumber ):
        Form( "Rootfiles/RawDataTiming-%08i-0000.root", runNumber );

    std::cout << "RawDataTiming - inputFile: " << inputFile << std::endl;
    std::cout << "RawDataTiming - pdfFile: " << pdfFile << std::endl;
    std::cout << "RawDataTiming - rootfileName: " << rootfileName << std::endl;

    // read hot channels
    if( use_mask )
    {
      const TString calibrationFile = Form( "DST/TPOT_HotChannels-%08i-0000.root", runNumber );
      std::cout << "RawDataTiming - calibrationFile: " << calibrationFile << std::endl;
      hot_channels.read( calibrationFile.Data() );
    }

    RootFile rootFile( rootfileName );
    PdfDocument pdfDocument( pdfFile );

    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "T" );

    // project adc vs sample vs strip in a 3D histogram
    // strips are grouped by chuncks of 64 corresponding to the 4 Resist region in each of the detectors
    const TString var( "samples.adc:samples.sample" );
    const TString var3d = Form( "%s:((samples.strip/64)+4*(samples.tile+8*(samples.layer-55)))", var.Data() );
    // auto h3d = new TH3I( "h3d", "h3d", 64, 0, 64, 150, 0, 150, 1200, 0, 1200 );
    auto h3d = new TH3I( "h3d", "h3d", 64, 0, 64, 360, 0, 360, 1200, 0, 1200 );
    h3d->GetXaxis()->SetTitle( "region_id" );
    h3d->GetYaxis()->SetTitle( "sample" );
    h3d->GetZaxis()->SetTitle( "adc" );

    if( use_mask )
    {
        Utils::max_entries = 200;
        const TCut hot_channel_cut = "is_good( samples.layer, samples.tile, samples.strip )";
        Utils::TreeToHisto( tree, h3d->GetName(), var3d, hot_channel_cut, false );
    } else {
        Utils::max_entries = 200;
        Utils::TreeToHisto( tree, h3d->GetName(), var3d, TCut(), false );
    }

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
            cv->Divide(2,2);

            for( int iregion = 0; iregion < 4; ++iregion )
            {
                // due to mapping internals, region are named opposite to strips.
                const int region = 4-iregion;
                const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
                const auto htitle = Form( "%s_R%i - %i,%i", name.c_str(), region, layer, itile );

                // get the bin matching layer, tile and region in the 3D histgoram
                const int bin = iregion + 4*(itile + 8*ilayer );
                h3d->GetXaxis()->SetRange( bin+1, bin+1 );

                // get the corresponding 2D histogram (adc vs sample)
                auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
                h2d->SetName( hname );
                h2d->SetTitle( htitle );

                // draw
                cv->cd(iregion+1);
                h2d->Draw("colz");
                rootFile.Add(h2d);
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
  { RawDataTiming( runnumber ); }
}
