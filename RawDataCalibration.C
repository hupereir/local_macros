#include <RootUtil/FileManager.h>

#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

#include <TProfile.h>
#include <map>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)


namespace
{
  // sample windows to determine signal and background clusters
  using sample_window_t = std::pair<unsigned short, unsigned short>;
  sample_window_t sample_window_background = {0, 20 };
}

TString RawDataCalibration( int runnumber = 26197 )
{

  // input
  const TString inputFile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataEvaluation-%08i-000?-full.root", runnumber );
  const TString calibrationFile = Form( "DST/TPOT_Pedestal-%08i-0000-test.root", runnumber );

  std::cout << "RawDataCalibration - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataCalibration - calibrationFile: " << calibrationFile << std::endl;
  
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  /// map fee id to Profile histogram
  using profile_map_t = std::map<unsigned short, TProfile*>;
  profile_map_t profile_map;

  // loop over tree entries
  const int entries = tree->GetEntries();
  // const int entries = 1000;
  std::cout << "RawDataClusterTree - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%100) )
    { std::cout << "RawDataCalibration - entry: " << i << std::endl; }

    tree->GetEntry(i);    

    // running profile map iterator
    profile_map_t::iterator piter = profile_map.begin();
    
    // loop over samples
    for( const auto& sample:container->samples )
    {
      // check timing
      if( sample.sample < sample_window_background.first || sample.sample >= sample_window_background.second )
      { continue; }

      const auto& fee = sample.fee_id;
      
      // see if fee id has changed since last check
      if( piter == profile_map.end() || piter->first != fee )
      {
      
        // find relevant profile histogram 
        piter= profile_map.lower_bound( fee );
        if( piter == profile_map.end() || fee < piter->first )
        {
          std::cout << "RawDataCalibration - creating profile for fee: " << fee << std::endl;
          
          // create and insert
          auto profile = new TProfile( Form( "h_adc_channel_%i", fee ), "ADC vs channel;channel;adc", MicromegasDefs::m_nchannels_fee, 0, MicromegasDefs::m_nchannels_fee );
          profile->SetErrorOption( "s" );
          piter = profile_map.insert(  piter, std::make_pair( fee, profile ) );      
        } 
      }
      
      // get corresponding profile and fill
      auto& profile = piter->second;
      profile->Fill( sample.channel, sample.adc );
    
    }
      
  }
  
  // create calibration data object
  MicromegasCalibrationData calibration_data;
  for( const auto& [fee, profile]:profile_map )
  {
    for( int i = 0; i < profile->GetNbinsX(); ++ i )
    {
      const auto pedestal = profile->GetBinContent(i+1);
      const auto rms = profile->GetBinError(i+1);
      calibration_data.set_pedestal( fee, i, pedestal );
      calibration_data.set_rms( fee, i, rms );
    }
  }
  calibration_data.write( calibrationFile.Data() );

  // print
  std::cout << calibration_data << std::endl;

  return calibrationFile;
}
