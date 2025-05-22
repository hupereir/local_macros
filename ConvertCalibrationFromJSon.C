R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasCalibrationData.h>
#include <nlohmann/json.hpp>

void ConvertCalibrationFromJSon( int runnumber = 45955 )
{

  const auto inputfile = Form("Calibrations/TPOT_Pedestal-%08i-0000.json", runnumber );
  const auto outputfile = Form("Calibrations/TPOT_Pedestal-%08i-0000-from_json.root", runnumber );

  std::cout << "EvaluateCalibration - inputfile: " << inputfile << std::endl;
  std::cout << "EvaluateCalibration - outputfile: " << outputfile << std::endl;

  std::ifstream in(inputfile);
  auto calib_list = nlohmann::json::parse(in);

  MicromegasCalibrationData m_calibration_data;

  static constexpr double pedestal = 60;
  static constexpr double nsigma = 5;

  for( const auto& calib:calib_list )
  {
    const auto fee_id = calib.at("fee_id").get<double>();
    const auto threshold = calib.at("threshold").get<double>();
    const double rms = (threshold/4-pedestal)/nsigma;
    std::cout
      << "fee_id: " << fee_id
      << " threshold: " << calib.at("threshold").get<int>()
      << " pedestal: " << pedestal
      << " rms: " << rms
      << std::endl;

    for( int ich = 0; ich < MicromegasDefs::m_nchannels_fee; ++ich )
    {
      // write to calibration data
      m_calibration_data.set_pedestal(fee_id, ich, pedestal );
      m_calibration_data.set_rms(fee_id, ich, rms );
    }

  }

  // write
  m_calibration_data.write(outputfile);

}
