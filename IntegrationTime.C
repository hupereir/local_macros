R__LOAD_LIBRARY(libRootUtilBase.so)
#include <RootUtil/Table.h>

// square
template<class T> T square( const T& x ) { return x*x; }

void IntegrationTime()
{

//   // rphi direction
//   constexpr std::array<double,10> segmentation = {{2e4, 1e4, 5e3, 2e3, 1e3, 500, 200, 100, 50, 0}};
//   constexpr double required = 150;
//
// //   // two layers
// //   constexpr std::array<double,10> residuals = {509, 586, 624, 665, 1005, 1574, 4378, 7714, 10800, 9000};
// //   constexpr double required = 150;
//
//   // one layer
//   constexpr std::array<double,9> residuals = {621, 628, 639, 704, 1407, 3005, 5830, 8238, 9549};

  // z direction
  constexpr std::array<double,10> segmentation = {{1.08e4, 5.4e3, 2.7e3, 1.08e3, 540, 270, 108, 54, 11, 0}};
  constexpr double required = 500;

//   // two layers
//   constexpr std::array<double,10> residuals = {{265, 266, 270, 433, 782, 1456, 2947, 3545, 3692, 3777 }};

  // one layer
  constexpr std::array<double,10> residuals = {{270, 271, 276, 556, 1043, 1907, 3427, 3721, 3633, 3777}};

  std::vector<double> timeArray;

  const double trk_frq = 400; // number of tracks per second per surface unit
  for( const auto residual:residuals )
  {
    const double nTrk = square( residual/required );
    const double time = nTrk/trk_frq;
    timeArray.push_back(time);
    std::cout << "IntegrationTime - residual: " << residual << " nTrk: " << nTrk << " integration: " << time << std::endl;
  }

  // print table

  Table table;
  table.AddColumn( "$N_\\rphi$", &segmentation[0], segmentation.size(), "%.0f" );
  table.AddColumn( "$\\sigma_\\rphi$ ($\\mu$m)", &residuals[0], residuals.size(), "%.0f" );
  table.AddColumn( "$\\Delta T$ (ms)", &timeArray[0], timeArray.size(), "%.3f" );
  table.ScaleLastColumn( 1000 );

  table.PrintLatex();

}
