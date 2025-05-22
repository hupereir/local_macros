#include <RootUtil/RootFile.h>
#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

// new bins (from Ross)
constexpr int m_phibins = 80;
constexpr int m_rbins = 52;

// divide by two because we store positive and negative z in different histograms
constexpr int m_zbins = 160 / 2;

// phi range
constexpr float m_phimin = 0;
constexpr float m_phimax = 2.*M_PI;

// r range
constexpr float m_rmin = 20;
constexpr float m_rmax = 78;

// z range
constexpr std::array<float,2> m_zmin = { -105.5, 0};
constexpr std::array<float,2> m_zmax = {0, 105.5};

// void CreateEmptyTimeOrderedDistortions()
// {
//   RootFile output( "distortion_maps/time_ordered_distortions_empty.root" );
// 
//   auto tree = new TTree( "TimeDists", "Tree of time-ordered distortions" );
// 
//   // create tree branches
//   int xingnum = 0;
//   tree->Branch( "xingnum", &xingnum );
//   
//   const std::array<TString, 3> hname =  
//   {
//     "hIntDistortionX",
//     "hIntDistortionY",
//     "hIntDistortionZ"
//   };
//   
//   const std::array<TString, 3> htitle = 
//   {
//     "Integrated X Distortion from (r,phi,z) to z=0 (centered in r,phi, and z)",
//     "Integrated Y Distortion from (r,phi,z) to z=0 (centered in r,phi, and z)",
//     "Integrated Z Distortion from (r,phi,z) to z=0 (centered in r,phi, and z)",
//   };
// 
//   for( int i = 0; i < 3; ++i )
//   {
//     auto h = new TH3F( hname[i], htitle[i],
//       m_phibins+2, m_phimin - (m_phimax-m_phimin)/m_phibins,  m_phimax + (m_phimax-m_phimin)/m_phibins, 
//       m_rbins+2, m_rmin - (m_rmax-m_rmin)/m_rbins,  m_rmax + (m_rmax-m_rmin)/m_rbins, 
//       m_zbins+2, m_zmin - (m_zmax-m_zmin)/m_zbins,  m_zmax + (m_zmax-m_zmin)/m_zbins );
//     h->GetXaxis()->SetTitle( "phi (rad)" );
//     h->GetYaxis()->SetTitle( "r (cm)" );
//     h->GetZaxis()->SetTitle( "z (cm)" );
//   
//     tree->Branch( hname[i], &h );
//   }
//   
//   const int nentries = 90;
//   for( int entry = 0; entry < nentries; ++entry )
//   { 
//     xingnum = entry;
//     tree->Fill(); 
//   }
//   
//   output.Add( tree );
//   
// }

void CreateEmptyStaticDistortions()
{
  
  RootFile output( "distortion_maps-new/empty.root" );
  
  const std::array<TString, 3> hname =  
  {
    "hIntDistortionP",
    "hIntDistortionR",
    "hIntDistortionZ"
  };
  
  std::array<TString, 2> suffix = { "_negz", "_posz" };
  
  
  for( int i = 0; i < 3; ++i )
  {
    for( int j =0; j<2; ++j )
    {
      
      const auto name = hname[i] + suffix[j];
      auto h = new TH3F( name, name,
        m_phibins+2, m_phimin - (m_phimax-m_phimin)/m_phibins,  m_phimax + (m_phimax-m_phimin)/m_phibins, 
        m_rbins+2, m_rmin - (m_rmax-m_rmin)/m_rbins,  m_rmax + (m_rmax-m_rmin)/m_rbins, 
        m_zbins+2, m_zmin[j] - (m_zmax[j]-m_zmin[j])/m_zbins,  m_zmax[j] + (m_zmax[j]-m_zmin[j])/m_zbins );
      h->GetXaxis()->SetTitle( "phi (rad)" );
      h->GetYaxis()->SetTitle( "r (cm)" );
      h->GetZaxis()->SetTitle( "z (cm)" );
      output.Add( h );
    }
  }
  
}


void CreateEmptyDistortions()
{
  CreateEmptyStaticDistortions();
  // CreateEmptyTimeOrderedDistortions();
}
