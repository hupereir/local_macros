#include <tpccalib/TpcSpaceChargeMatrixContainerv1.h>
#include <tpccalib/TpcSpaceChargeMatrixInversion.h>
#include <RootUtil/FileManager.h>
#include <g4eval_hp/TrackingEvaluator_hp.h>

#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>

#include <Eigen/Dense>

R__LOAD_LIBRARY(libg4eval_hp.so)
R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libtpccalib.so)

#include "LayerDefines.h"

namespace
{

  template< class T>
    constexpr T square( const T& x ) { return x*x; }

  template< class T>
    T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  // matrix container for distortion correction
  TpcSpaceChargeMatrixContainerv1 matrix_container;

}


// phi range
constexpr float m_phimin = 0;
constexpr float m_phimax = 2.*M_PI;

// r range
constexpr float m_rmin = 20;
constexpr float m_rmax = 78;

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;

static constexpr bool use_truth_information = false;
static constexpr bool use_micromegas = true;

// histogram limits
static constexpr float max_residual_drphi = 0.5;
static constexpr float max_residual_dz = 0.5;

// static constexpr float max_residual_drphi = 1.0;
// static constexpr float max_residual_dz = 1.0;

// // larger values to handle full distortions
// static constexpr float max_residual_drphi = 2.5;
// static constexpr float max_residual_dz = 2.5;

static constexpr float max_alpha = 0.5;
static constexpr float max_beta = 1.0;

static constexpr int min_entries = 100;

//__________________________________________________________________
bool accept_track( const TrackingEvaluator_hp::TrackStruct& track )
{

  // momentum cut
  if( use_truth_information )
  {

    if( track._truth_pt < 0.5 ) return false;

  } else {

    if( track._pt < 0.5 ) return false;

  }

  // hit pattern cuts
  if( track._nclusters_mvtx < 3 ) return false;
  if( track._nclusters_intt < 2 ) return false;
  if( use_micromegas && track._nclusters_micromegas < 2 ) return false;

//   if( use_micromegas )
//   {
//     // loop over clusters
//     // restrict to tiles 0-3
//     for( const auto& cluster:track._clusters )
//     {
//       if( (cluster._layer == 55 || cluster._layer == 56) && cluster._tileid >= 4 )
//       { return false; }
//     }
//   }

  // cut on charge
  // if( !(track._charge>0) ) return false;
  // if( !(track._charge<0) ) return false;

  return true;
}

//__________________________________________________________________
bool accept_cluster( const TrackingEvaluator_hp::ClusterStruct& cluster )
{

  // skip clusters not in TPC
  if( cluster._layer < firstLayer_tpc || cluster._layer >= firstLayer_tpc + nLayers_tpc )
  { return false; }

  if( true )
  {
    /*
    remove clusters with too small errors since they are likely pathological
    and have a large contribution to the chisquare
    */
    if( cluster._trk_r*cluster._phi_error < 0.015 ) return false;
    if( cluster._z_error < 0.05 ) return false;
  }

  // check against window
  if( use_truth_information )
  {
    if( std::abs( cluster._truth_alpha ) > max_alpha ) return false;
    if( std::abs( cluster._truth_beta ) > max_beta ) return false;
  } else {
    if( std::abs( cluster._trk_alpha ) > max_alpha ) return false;
    if( std::abs( cluster._trk_beta ) > max_beta ) return false;
  }

  return true;
}

//_________________________________________________________________________
int get_cell( const TrackingEvaluator_hp::ClusterStruct& cluster )
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  matrix_container.get_grid_dimensions( phibins, rbins, zbins );

  // phi
  float phi = cluster._phi;
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  const int iphi = phibins * (phi - m_phimin) / (m_phimax - m_phimin);

  // r
  const auto r = cluster._r;
  if( r < m_rmin || r >= m_rmax ) return -1;
  const int ir = rbins * (r - m_rmin) / (m_rmax - m_rmin);

  // z
  const auto z = cluster._z;
  if( z < m_zmin || z >= m_zmax ) return -1;
  const int iz = zbins * (z - m_zmin) / (m_zmax - m_zmin);

  // get index from matrix container
  return matrix_container.get_cell_index( iphi, ir, iz );
}

//_________________________________________________________________________
void Distortions_parallel(
  const char* inputfile,
  const char* outputfilename
  )
{

  std::cout << "Distortions - inputfile: " << inputfile << std::endl;
  std::cout << "Distortions - outputfilename: " << outputfilename << std::endl;

  // reconstructed distortion grid size
  // matrix_container.set_grid_dimensions( 36, 16, 80 );
  // matrix_container.set_grid_dimensions( 36, 48, 80 );
  matrix_container.set_grid_dimensions( 180, 16, 80 );

  // get tree
  FileManager fileManager( inputfile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return;

  auto container = new TrackingEvaluator_hp::Container;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  // counters
  int total_tracks = 0;
  int accepted_tracks = 0;

  int total_clusters = 0;
  int accepted_clusters = 0;

  // loop over entries
  const int entries = tree->GetEntries();
  const int offset = 0;
  std::cout << "Distortions - entries: " << entries << " (" << tree->GetEntries() << ")" << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    tree->GetEntry(i + offset);
    if( !(i%100) )
    { std::cout << "Distortions - entry: " << i << std::endl;}

    // loop over tracks
    for( const auto& track:container->tracks() )
    {

      ++total_tracks;

      // check
      if( !accept_track( track ) ) continue;
      ++accepted_tracks;

      // loop over clusters
      for( const auto& cluster:track._clusters )
      {

        ++total_clusters;
        if( !accept_cluster( cluster ) ) continue;

        // store errors
        const auto erp = use_truth_information ?
          square(cluster._trk_r*cluster._phi_error):
          ( square(cluster._trk_r)*( square( cluster._phi_error ) + square( cluster._trk_phi_error ) ) );

        const auto ez = use_truth_information ?
          square(cluster._z_error):
          ( square( cluster._z_error ) + square( cluster._trk_z_error ) );

        // sanity check
        if( std::isnan( erp ) ) continue;
        if( std::isnan( ez ) ) continue;

        const auto r = use_truth_information ?  cluster._truth_r:cluster._trk_r;
        const auto drp = r * (use_truth_information ? delta_phi( cluster._phi - cluster._truth_phi ):delta_phi( cluster._phi - cluster._trk_phi ));
        if( std::isnan(drp) ) continue;

        const auto talpha = use_truth_information ?
          -std::tan( cluster._truth_alpha ):
          -std::tan( cluster._trk_alpha );
        if( std::isnan(talpha) ) continue;

        const auto dz = use_truth_information ? (cluster._z - cluster._truth_z):(cluster._z - cluster._trk_z );
        if( std::isnan(dz) ) continue;

        const auto tbeta = use_truth_information ?
          -std::tan( cluster._truth_beta ):
          -std::tan( cluster._trk_beta );
        if( std::isnan(tbeta) ) continue;

        // check against max residuals
        if( std::abs( drp ) > max_residual_drphi ) continue;
        if( std::abs( dz ) > max_residual_dz ) continue;

        // find relevant cell
        auto index = get_cell( cluster );
        if( index < 0 ) continue;

        // Fill distortion matrices
        matrix_container.add_to_lhs(index, 0, 0, 1./erp );
        matrix_container.add_to_lhs(index, 0, 1, 0 );
        matrix_container.add_to_lhs(index, 0, 2, talpha/erp );

        matrix_container.add_to_lhs(index, 1, 0, 0 );
        matrix_container.add_to_lhs(index, 1, 1, 1./ez );
        matrix_container.add_to_lhs(index, 1, 2, tbeta/ez );

        matrix_container.add_to_lhs(index, 2, 0, talpha/erp );
        matrix_container.add_to_lhs(index, 2, 1, tbeta/ez );
        matrix_container.add_to_lhs(index, 2, 2, square(talpha)/erp + square(tbeta)/ez );

        matrix_container.add_to_rhs(index, 0, drp/erp );
        matrix_container.add_to_rhs(index, 1, dz/ez );
        matrix_container.add_to_rhs(index, 2, talpha*drp/erp + tbeta*dz/ez );

        // update entries in cell
        matrix_container.add_to_entries(index);

        // increment number of accepted clusters
        ++accepted_clusters;

      }
    }
  }

  // print statistics
  std::cout << "Distortions - entries: " << entries << std::endl;
  std::cout << "Distortions - track statistics total: " << total_tracks << " accepted: " << accepted_tracks << " fraction: " << 100.*accepted_tracks/total_tracks << "%" << std::endl;
  std::cout << "Distortions - cluster statistics total: " << total_clusters << " accepted: " << accepted_clusters << " fraction: " << 100.*accepted_clusters/total_clusters << "%" << std::endl;

  {
    std::unique_ptr<TFile> outputfile( TFile::Open( outputfilename, "RECREATE" ) );
    outputfile->cd();
    matrix_container.Write( "TpcSpaceChargeMatrixContainer" );
  }

  return;
}

//_________________________________________________________________________
void MatrixInversion( const char* inputfile, const char* outputfilename )
{
  // do the insversion and write to file
  TpcSpaceChargeMatrixInversion matrix_inversion;
  matrix_inversion.set_outputfile( outputfilename );

  FileManager f( inputfile );
  for( const auto& filename:f.GetFiles() )
  { matrix_inversion.add_from_file( filename.Data()); }

  matrix_inversion.calculate_distortions();
}

//_________________________________________________________________________
void Distortions_parallel()
{
  // const TString tag = "_flat_acts_truth_no_distortion";
  const TString tag = "_flat_acts_truth_distorted";
  const TString inputfile = Form( "DST/CONDOR%s/dst_reco%s_?.root", tag.Data(), tag.Data() );

  const TString matrix_container_file = Form( "DST/CONDOR%s/SpaceChargeMatrixContainer%s.root", tag.Data(), tag.Data() );
  Distortions_parallel( inputfile, matrix_container_file );

  // root file
  TString outputfile;
  if( use_truth_information ) outputfile = Form( "Rootfiles/Distortions_full%s_truth.root", tag.Data() );
  else if( use_micromegas ) outputfile = Form( "Rootfiles/Distortions_full%s_mm.root", tag.Data() );
  else outputfile = Form( "Rootfiles/Distortions_full%s_all.root", tag.Data() );
  MatrixInversion( matrix_container_file, outputfile );
}

//_________________________________________________________________________
void MatrixInversion( TString tag = TString() )
{
  if( tag.IsNull() )
  {
    // tag = "_flat_acts_truth_nodistortion";
    // tag = "_flat_acts_truth_distorted";
    // tag = "_flat_acts_truth_notpc_nodistortion";
    // tag = "_flat_acts_truth_notpc_distorted";
    // tag = "_flat_genfit_truth_notpc_nodistortion-new";
    // tag = "_flat_genfit_truth_notpc_nodistortion-nosurvey-new";
    tag = "_flat_genfit_truth_notpc_nodistortion-fix2";
    // tag = "_flat_genfit_truth_notpc_distorted";
    // tag = "_flat_genfit_truth_notpc_distorted-nosurvey";
    // tag = "_flat_genfit_truth_notpc_distorted-fix";
  }

  // const TString inputfile = Form( "DST/CONDOR%s/SpaceChargeMatrixContainer%s_*.root", tag.Data(), tag.Data() );
  const TString inputfile = Form( "DST/CONDOR%s/TpcSpaceChargeMatrices_offline_test_%s_*.root", tag.Data(), tag.Data() );

  // root file
  TString outputfile;
  if( use_truth_information ) outputfile = Form( "Rootfiles/Distortions_offline_full%s_truth.root", tag.Data() );
  else if( use_micromegas ) outputfile = Form( "Rootfiles/Distortions_offline_full%s_mm.root", tag.Data() );
  else outputfile = Form( "Rootfiles/Distortions_full%s_all.root", tag.Data() );
  MatrixInversion( inputfile, outputfile );
}

