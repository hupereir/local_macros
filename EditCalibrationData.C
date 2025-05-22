#include <micromegas/MicromegasCalibrationData.h>

R__LOAD_LIBRARY(libmicromegas.so)

void EditCalibrationData()
{

  MicromegasCalibrationData calibData;
  // calibData.read( "DST/TPOT_Pedestal-00009416-0000.root" );
  calibData.read( "DST/TPOT_Pedestal-00026174-0000.root" );
  
  if( true )
  {
    // fee, channel, pedestal, rms
    using calibration_data_t = std::tuple<int, int, double, double>;
    using calibration_vector_t = std::vector<calibration_data_t>;
    calibration_vector_t calibrations = {
      { 12, 165, 147.1, 5.935 },
      { 18, 110, 70.2, 7.4 },
      { 18, 117, 75.6, 7.3 },
      { 18, 127, 87.7, 7.63 },
      { 18, 139, 76.8, 7.59 },
      { 18, 139, 75.9, 7.66 }
    };
    
    for( const auto& [fee, channel, pedestal, rms]:calibrations )
    {    
      calibData.set_pedestal( fee, channel, pedestal );
      calibData.set_rms( fee, channel, rms );
    }
  }
  calibData.write( "DST/TPOT_Pedestal-00009416-0000.root" );
  
  std::cout << calibData << std::endl;
  
}
