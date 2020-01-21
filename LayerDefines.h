#ifndef LayerDefines_h
#define LayerDefines_h

#include <array>

// convenient defines
constexpr int nDetectors = 6;
constexpr int nLayersTotal = 57;

constexpr std::array<int, nDetectors> nLayers = {{ 3, 4, 16, 16, 16, 2 }};
constexpr std::array<int, nDetectors> firstLayer = {{ 0, 3, 7, 23, 39, 55 }};
constexpr std::array<float,nLayersTotal> radius = {{
  2.41, 3.17, 4.05, 9.01, 9.45, 11.25, 11.25, 30.15, 31.05, 31.95,
  31.95, 32.85, 33.75, 33.75, 34.65, 35.55, 35.55, 36.45, 37.35, 38.25,
  38.25, 39.15, 40.05, 40.95, 41.85, 42.75, 44.55, 45.45, 47.25, 48.15,
  49.05, 50.85, 51.75, 53.55, 54.45, 55.35, 57.15, 58.05, 58.95, 60.75,
  61.65, 62.55, 63.45, 64.35, 66.15, 67.05, 67.95, 68.85, 69.75, 71.55,
  72.45, 73.35, 74.25, 75.15, 76.05, 82.35, 84.15
}};

int get_detid( int layer )
{
  for( int idet = 0; idet < nDetectors; ++idet )
  { if( layer < firstLayer[idet] ) return idet-1; }

  return layer < nLayersTotal ? nDetectors-1:nDetectors;
}

#endif
