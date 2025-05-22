#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>

R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)

void GetGeometry()
{
  // reco const
  auto rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", 6);
  const auto geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  std::cout << "GetGeometry - geofile: " << geofile << std::endl;

  const auto alignmentParamsFile = CDBInterface::instance()->getUrl("TRACKINGALIGNMENT");
  std::cout << "GetGeometry - alignmentParamsFile: " << alignmentParamsFile << std::endl;
}
