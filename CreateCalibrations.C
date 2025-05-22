#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>

R__LOAD_LIBRARY(libmicromegas.so)


void CreateCalibrations()
{

  const std::string calibration_filename = "TPOT_Pedestal-000-test.root";
  std::cout << "CreateCalibrations - calibration_filename: " << calibration_filename << std::endl;

  /*
   * define default pedestal and default RMS
   * threshold is set to pedestal + 5x rms
   */
  static constexpr double pedestal = 70;
  static constexpr double rms = 10;

  // create calibration data object
  MicromegasCalibrationData calibration_data;

  // micromegas mapping
  MicromegasMapping mapping;

  // loop over fees
  for( const auto& fee_id:mapping.get_fee_id_list() )
  {

    // loop over channels
    for( int i = 0; i < MicromegasDefs::m_nchannels_fee; ++i )
    {
      calibration_data.set_pedestal( fee_id, i, pedestal );
      calibration_data.set_rms( fee_id, i, rms );
    }
  }

  calibration_data.write( calibration_filename );

}
