#ifndef LayerDefines_h
#define LayerDefines_h

#include <array>

// convenient defines
constexpr int nDetectors = 7;
constexpr int nLayersTotal = 57;

inline constexpr bool is_tpc( int idet ) { return idet >= 2 && idet < 5; }

constexpr int nLayers_tpc = 48;
constexpr int firstLayer_tpc = 7;

constexpr float rmin_tpc = 20;
constexpr float rmax_tpc = 78;
constexpr float rlength_tpc = rmax_tpc - rmin_tpc;

constexpr std::array<int, nDetectors> nLayers = {{ 3, 4, 16, 16, 16, 1, 1 }};
constexpr std::array<int, nDetectors> firstLayer = {{ 0, 3, 7, 23, 39, 55, 56 }};
constexpr std::array<double,nLayersTotal> radius = {{
  2.58, 3.37, 4.15,
  9.01, 9.56, 10.85, 11.38,
  30.32, 30.93, 31.57, 32.19, 32.81, 33.44, 34.06, 34.69, 35.31, 35.94, 36.56, 37.19, 37.81, 38.44, 39.06, 39.69,
  40.62, 41.88, 43.12, 44.38, 45.62, 46.88, 48.12, 49.38, 50.62, 51.88, 53.12, 54.38, 55.62, 56.88, 58.12, 59.38,
  60.53, 61.60, 62.65, 63.72, 64.79, 65.85, 66.90, 67.96, 69.04, 70.10, 71.15, 72.21, 73.29, 74.35, 75.40, 76.46,
  82.00, 84.00}};

constexpr std::array<double,nLayers_tpc+1> tpc_radius = 
{{
  30, 30.625, 31.25, 31.875, 32.5, 33.125, 33.75, 34.375, 35, 35.625, 
  36.25, 36.875, 37.5, 38.125, 38.75, 39.375, 40, 41.25, 42.5, 43.75, 
  45, 46.25, 47.5, 48.75, 50, 51.25, 52.5, 53.75, 55, 56.25, 
  57.5, 58.75, 60, 61.0625, 62.125, 63.1875, 64.25, 65.3125, 66.375, 67.4375, 
  68.5, 69.5625, 70.625, 71.6875, 72.75, 73.8125, 74.875, 75.9375, 77 
}};
    
int get_detid( int layer )
{
  for( int idet = 0; idet < nDetectors; ++idet )
  { if( layer < firstLayer[idet] ) return idet-1; }

  return layer < nLayersTotal ? nDetectors-1:nDetectors;
}

#endif
