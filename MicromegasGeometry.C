R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <micromegas/MicromegasMapping.h>
#include <g4eval_hp/MicromegasGeometryContainer.h>

namespace
{

  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  template<class T>
    inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  // geometry container
  std::unique_ptr<MicromegasGeometryContainer> m_geometry_container;

  // mapping
  MicromegasMapping mapping;

  // detector names
  std::vector<std::string> detector_names;

  // load geometry from file
  void load_geometry( const std::string& filename )
  {
    auto inputfile( TFile::Open( filename.c_str() ) );
    m_geometry_container.reset( dynamic_cast<MicromegasGeometryContainer*>( inputfile->Get( "MicromegasGeometryContainer" ) ) );
    m_geometry_container->identify();
  }

  // save detector names
  void save_detector_names()
  {

    // get detector names that match tile/layer
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name ) );
    }
  }

  template<class T, class U>
    std::ostream& operator << (std::ostream& out, const std::pair<T, U> pair)
  {
    out << "{" << pair.first << "," << pair.second << "}";
    return out;
  }

  // generic container printout
  template<template<typename,typename> class Container, class T, class A>
    std::ostream& operator << (std::ostream& out, const Container<T,A>& container)
  {
    out << "{";
    bool first = true;
    for( const auto& value:container )
    {
      if( !first ) { out << ", "; };
      first = false;
      out << value;
    }
    out << "}";
    return out;
  }

}

void MicromegasGeometry()
{

  // load geometry
  load_geometry( "micromegas_geometry.root" );

  // save detector names
  save_detector_names();

  if( true )
  {
    // get minimum radius of each chamber
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      //  int ilayer = 0;
      for( int tile = 0; tile < 8; ++tile )
    {
      int detid = tile + 8*(ilayer);
      int layer = 55+ilayer;

      // loop over strips
      double r_min = -1;
      int strip_min = -1;
      for( int strip=0; strip < 256; ++strip )
      {
        const auto strip_begin = m_geometry_container->get_strip_begin( layer, tile, strip );
        const auto strip_end = m_geometry_container->get_strip_end( layer, tile, strip );

        const auto strip_center = (strip_begin + strip_end)*0.5;
        const auto radius = get_r( strip_center.x(), strip_center.y() );

        if( r_min <0 || radius < r_min )
        {
          r_min = radius;
          strip_min = strip;
        }
      }

      std::cout << "detector: " << detector_names[detid] << " layer: " << layer << " tile: " << tile << " r_min: " << r_min << " strip_min: " << strip_min << std::endl;
    }
  }

  using phi_window_t = std::pair<double, double>;
  using phi_window_list_t = std::vector<phi_window_t>;
  phi_window_list_t phi_window_list(8,phi_window_t{});

  if( true )
  {

    // get phi range of each module
    const int ilayer = 0;
    for( int tile = 0; tile < 8; ++tile )
    {
      int detid = tile + 8*(ilayer);
      int layer = 55+ilayer;

      auto get_strip_phi = [layer, tile]( int strip )
      {
        const auto strip_begin = m_geometry_container->get_strip_begin( layer, tile, strip );
        const auto strip_end = m_geometry_container->get_strip_end( layer, tile, strip );
        const auto strip_center = (strip_begin + strip_end)*0.5;
        return std::atan2(strip_center.y(), strip_center.x());
      };

      phi_window_t phi_window( get_strip_phi(0), get_strip_phi(255) );
      phi_window_list[tile] = phi_window;
      // std::cout << "detector: " << detector_names[detid] << " layer: " << layer << " tile: " << tile << " phi_window=" << phi_window << std::endl;
    }

    auto get_phi_window = [phi_window_list]( std::vector<int> indexes )
    {
      phi_window_t out = {0,0};
      for(const auto& i:indexes)
      {
        const auto& window = phi_window_list[i];
        if( out.first == 0 || window.first > out.first ) out.first = window.first;
        if( out.second == 0 || window.second < out.second ) out.second = window.second;
      }
      return out;
    };

    const auto phi_window_central = get_phi_window( {0,1,2,3} );
    const auto phi_window_east = get_phi_window( {4,5} );
    const auto phi_window_west = get_phi_window( {6,7} );
    std::cout << "phi_window_central=" << phi_window_central << std::endl;
    std::cout << "phi_window_east=" << phi_window_east << std::endl;
    std::cout << "phi_window_west=" << phi_window_west << std::endl;

  }

  using theta_window_t = std::pair<double, double>;
  using theta_window_list_t = std::vector<theta_window_t>;
  theta_window_list_t theta_window_list(8, theta_window_t{});

  using eta_window_t = std::pair<double, double>;
  using eta_window_list_t = std::vector<eta_window_t>;
  eta_window_list_t eta_window_list(8, eta_window_t{});

  if( true )
  {

    // get phi range of each module
    const int ilayer = 1;
    for( int tile = 0; tile < 8; ++tile )
    {
      int detid = tile + 8*(ilayer);
      int layer = 55+ilayer;

      auto get_strip_theta = [layer, tile]( int strip )
      {
        const auto strip_begin = m_geometry_container->get_strip_begin( layer, tile, strip );
        const auto strip_end = m_geometry_container->get_strip_end( layer, tile, strip );
        const auto strip_center = (strip_begin + strip_end)*0.5;
        // return std::atan2(strip_center.z(), get_r(strip_center.x(), strip_center.y()));
        return std::atan2(get_r(strip_center.x(), strip_center.y()), strip_center.z());
      };

      auto get_eta = []( double theta ) { return -std::log(std::tan(theta/2)); };

      const theta_window_t theta_window(get_strip_theta(0),get_strip_theta(255));
      const eta_window_t eta_window{ get_eta( theta_window.first ), get_eta( theta_window.second ) };
      theta_window_list[tile] = theta_window;
      eta_window_list[tile] = eta_window;
      // std::cout << "detector: " << detector_names[detid] << " layer: " << layer << " tile: " << tile << " theta_window=" << theta_window << std::endl;
    }

    std::cout << "theta_window_central=" << theta_window_list_t{theta_window_list[0], theta_window_list[1], theta_window_list[2], theta_window_list[3]} << std::endl;
    std::cout << "theta_window_east=" << theta_window_list_t{theta_window_list[4], theta_window_list[5]} << std::endl;
    std::cout << "theta_window_west=" << theta_window_list_t{theta_window_list[6], theta_window_list[7]} << std::endl;
  }

  {
    // print all tile information
    const int ilayer = 0;
    const int layer = ilayer+55;
    for( int itile = 0; itile < 8; ++itile )
    {
      int detid = itile+8*ilayer;
      std::cout << "detector: " << detector_names[detid] << " layer: " << layer << " tile: " << itile
        << " phi_window=" << phi_window_list[itile]
        << " (rad) theta_window=" << theta_window_list[itile]
        << " (rad) eta_window=" << eta_window_list[itile]
        << std::endl;
    }
  }

}
