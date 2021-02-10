#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

#include <Eigen/Dense>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

// binning and range
constexpr int m_phibins = 36;
constexpr int m_rbins = 16;
constexpr int m_zbins = 80;

// // new bins (from Ross)
// constexpr int m_phibins = 80;
// constexpr int m_rbins = 52;
// constexpr int m_zbins = 160;

// phi range
constexpr float m_phimin = 0;
constexpr float m_phimax = 2.*M_PI;

// r range
constexpr float m_rmin = 20;
constexpr float m_rmax = 78;

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;

constexpr int m_totalbins = m_zbins*m_phibins*m_rbins;

static constexpr bool use_truth_information = false;
static constexpr bool use_micromegas = true;
static constexpr bool save_histograms = false;

// specify bins for which 2D histograms will be saved
// static constexpr int phibin_rec = 11;
// static constexpr int zbin_rec = 42;

static constexpr int phibin_rec = 17;
static constexpr int zbin_rec = 53;

// histogram limits
// static constexpr float max_residual_drphi = 0.5;
// static constexpr float max_residual_dz = 0.5;

// larger values to handle full distortions
static constexpr float max_residual_drphi = 2.0;
static constexpr float max_residual_dz = 2.0;

static constexpr float max_alpha = 0.5;
static constexpr float max_beta = 1.0;

//_________________________________________________________________________
template< class T>
constexpr T square( const T& x ) { return x*x; }

//_________________________________________________________________________
template< class T>
T delta_phi( const T& phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//__________________________________________________________________
bool accept_track( const TrackingEvaluator_hp::TrackStruct& track )
{
  
  // momentum cut
  if( track._truth_pt < 0.5 ) return false;
  if( track._pt < 0.5 ) return false;
  
  // hit pattern cuts
  if( track._nclusters_mvtx < 3 ) return false;
  if( track._nclusters_intt < 2 ) return false;
  if( use_micromegas && track._nclusters_micromegas < 2 ) return false;
  
  // also cut on the ndf
  // if( use_micromegas && track._ndf < 9 ) return false;

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
int get_cell( int iphi, int ir, int iz )
{
  if( iphi < 0 || iphi >= m_phibins ) return -1;
  if( ir < 0 || ir >= m_rbins ) return -1;
  if( iz < 0 || iz >= m_zbins ) return -1;
  return iz + m_zbins*( ir + m_rbins*iphi );
}

//_________________________________________________________________________
int get_cell( float phi, float r, float z )
{

  // phi
  // bound check
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  int iphi = m_phibins*(phi-m_phimin)/(m_phimax-m_phimin);

  // radius
  if( r < m_rmin || r >= m_rmax ) return -1;
  int ir = m_rbins*(r-m_rmin)/(m_rmax-m_rmin);

  // z
  if( z < m_zmin || z >= m_zmax ) return -1;
  int iz = m_zbins*(z-m_zmin)/(m_zmax-m_zmin);

  return get_cell( iphi, ir, iz );
}

//_________________________________________________________________________
int get_cell( const TrackingEvaluator_hp::ClusterStruct& cluster )
{ return get_cell( cluster._phi, cluster._r, cluster._z ); }

//_________________________________________________________________________
double get_correction( TH3* h, const TrackingEvaluator_hp::ClusterStruct& cluster )
{
  
  // phi range
  constexpr float m_phimin = 0;
  constexpr float m_phimax = 2.*M_PI;

  // phi
  auto phi = cluster._phi;
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  
  auto iphi = h->GetXaxis()->FindBin( phi );
  if( iphi == 0 || iphi > h->GetNbinsX() ) return 0;
  
  // radius
  auto ir = h->GetYaxis()->FindBin( cluster._r );
  if( ir == 0 || ir > h->GetNbinsY() ) return 0;
  
  // z
  auto iz = h->GetZaxis()->FindBin( cluster._z );
  if( iz == 0 || iz > h->GetNbinsZ() ) return 0;
  
  return h->GetBinContent( iphi, ir, iz );
}

//_________________________________________________________________________
TString Distortions_iterate()
{
  std::cout << "Distortions_iterate - bins: " << m_totalbins << std::endl;

  // warning: needs update when using map with larger distortions
  using TH1_ptr_t = std::unique_ptr<TH1>;
  using TH2_ptr_t = std::unique_ptr<TH2>;

  // histogram arrays
  std::vector<TH1_ptr_t> residuals_drphi( m_totalbins );
  std::vector<TH1_ptr_t> residuals_dz( m_totalbins );
  std::vector<TH2_ptr_t> residuals_2d_drphi( m_totalbins );
  std::vector<TH2_ptr_t> residuals_2d_dz( m_totalbins );

  if( save_histograms )
  {
    for( int iphi = 0; iphi < m_phibins; ++iphi )
      for( int ir = 0; ir < m_rbins; ++ir )
      for( int iz = 0; iz < m_zbins; ++iz )
    {
      const auto icell = get_cell( iphi, ir, iz );

      {
        // rphi residuals
        const auto hname = Form( "residual_drphi_p%i_r%i_z%i", iphi, ir, iz );
        residuals_drphi[icell].reset( new TH1F( hname, hname, 100, -max_residual_drphi, +max_residual_drphi ) );
        residuals_drphi[icell] ->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
      }

      if( iphi+1 == phibin_rec && iz+1 == zbin_rec )
      {
        // 2D histograms
        const auto hname = Form( "residual_2d_drphi_p%i_r%i_z%i", iphi, ir, iz );
        residuals_2d_drphi[icell].reset( new TH2F( hname, hname, 100, -std::tan(max_alpha), std::tan(max_alpha), 100, -max_residual_drphi, +max_residual_drphi ) );
        residuals_2d_drphi[icell]->GetXaxis()->SetTitle( "tan#alpha" );
        residuals_2d_drphi[icell]->GetYaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
      }

      {
        // z residuals
        const auto hname = Form( "residual_dz_p%i_r%i_z%i", iphi, ir, iz );
        residuals_dz[icell].reset( new TH1F( hname, hname, 100, -max_residual_dz, +max_residual_dz ) );
        residuals_dz[icell]->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
      }

      if( iphi+1 == phibin_rec && iz+1 == zbin_rec )
      {
        // 2D histograms
        const auto hname = Form( "residual_2d_dz_p%i_r%i_z%i", iphi, ir, iz );
        residuals_2d_dz[icell].reset( new TH2F( hname, hname, 100, -std::tan(max_beta), std::tan(max_beta), 100, -max_residual_dz, +max_residual_dz ) );
        residuals_2d_dz[icell]->GetXaxis()->SetTitle( "tan#alpha" );
        residuals_2d_dz[icell]->GetYaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
      }

    }
  }

  // input files
  const TString tag = "_realistic_micromegas";
  const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc_distortions/dst_reco_realistic_micromegas_*.root";

  const TString subtag = "-test-loose";
  
  // root file
  TString rootfilename;
  if( use_truth_information ) rootfilename = Form( "Rootfiles/Distortions_iterate%s_truth%s.root", tag.Data(), subtag.Data() );
  else if( use_micromegas ) rootfilename = Form( "Rootfiles/Distortions_iterate%s_mm%s.root", tag.Data(), subtag.Data() );
  else rootfilename = Form( "Rootfiles/Distortions_iterate%s_all%s.root", tag.Data(), subtag.Data() );
  RootFile rootFile( rootfilename );

  // log file
  TString logfilename;
  if( use_truth_information ) logfilename = Form( "Distortions_iterate%s_truth%s.log", tag.Data(), subtag.Data() );
  else if( use_micromegas ) logfilename = Form( "Distortions_iterate%s_mm%s.log", tag.Data(), subtag.Data() );
  else logfilename = Form( "Distortions_iterate%s_all%s.log", tag.Data(), subtag.Data() );
  std::ofstream out( logfilename.Data() );

  std::cout << "Distortions_iterate - inputfile: " << inputFile << std::endl;
  std::cout << "Distortions_iterate - rootfile: " << rootfilename << std::endl;
  std::cout << "Distortions_iterate - logfile: " << logfilename << std::endl;

  // get tree
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return rootfilename;

  auto container = new TrackingEvaluator_hp::Container;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  // distortion correction histograms
  auto hentries_rec = new TH3F( "hentries_rec", "hentries_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hphi_rec = new TH3F( "hDistortionP_rec", "hDistortionP_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hr_rec = new TH3F( "hDistortionR_rec", "hDistortionR_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hz_rec = new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );

  // set axis labels
  for( auto h:{ hentries_rec, hphi_rec, hr_rec, hz_rec } )
  {
    h->GetXaxis()->SetTitle( "#phi (rad)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "z (cm)" );
    rootFile.Add( h );
  }
  
  // define iterations
  static constexpr int niterations = 3;
  static constexpr std::array<double,niterations> nsigmacut = {{0, 3, 2}};
  
  // pull histograms
  std::vector<TH1_ptr_t> pulls_drphi( niterations );
  std::vector<TH1_ptr_t> pulls_dz( niterations );
  for( int iiter = 0; iiter < niterations; ++iiter )
  {
    {
      const auto hname = Form( "pulls_drphi_%i", iiter );
      pulls_drphi[iiter].reset( new TH1F( hname, hname, 100, -10, 10 ) );
      pulls_drphi[iiter]->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} pulls" );
    }      
    
    {
      const auto hname = Form( "pulls_dz_%i", iiter );
      pulls_dz[iiter].reset( new TH1F( hname, hname, 100, -10, 10 ) );
      pulls_dz[iiter]->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} pulls" );
    }      
  }
  
  // counters
  int total_tracks = 0;  
  int total_clusters = 0;
  int accepted_tracks = 0;
  std::array<int,niterations> accepted_clusters;
  
  // loop over entries
  // const int entries = tree->GetEntries();
  const int entries = 500000;
  // const int entries = 4500000;
  std::cout << "Distortions_iterate - entries: " << entries << std::endl;
  
  // iteration loop
  for( int iiter = 0; iiter < niterations; ++iiter )
  {

    // reset all histograms
    for( const auto& h:residuals_drphi ) { if( h ) h->Reset(); }
    for( const auto& h:residuals_dz ) { if( h ) h->Reset(); }
    for( const auto& h:residuals_2d_drphi ) { if( h ) h->Reset(); }
    for( const auto& h:residuals_2d_dz ) { if( h ) h->Reset(); }
      
    // reset counter
    total_tracks = 0;
    total_clusters = 0;

    accepted_tracks = 0;
    accepted_clusters[iiter] = 0;

    // entries and matrices
    // state vector is s = (deltarphi, deltaz, deltar0) (3 rows, 1 column)
    // equation is lhs.s = rhs
    static constexpr int ncoord = 3;
    using matrix_t = Eigen::Matrix<float, ncoord, ncoord >;
    using column_t = Eigen::Matrix<float, ncoord, 1 >;
        
    std::vector<int> binentries( m_totalbins, 0 );
    std::vector<matrix_t> lhs( m_totalbins, matrix_t::Zero() );
    std::vector<column_t> rhs( m_totalbins, column_t::Zero() );
    
    for( int i = 0; i < entries; ++i )
    {
      tree->GetEntry(i);
      if( !(i%100) )
      { std::cout << "Distortions_iterate - iteration: " << iiter << " entry: " << i << std::endl; }
      
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

          // rphi residual
          const auto drp = use_truth_information ?
            cluster._truth_r*( delta_phi( cluster._phi - cluster._truth_phi ) ):
            cluster._trk_r*( delta_phi( cluster._phi - cluster._trk_phi ) );
          if( std::isnan(drp) ) continue;

          // rphi angle
          const auto talpha = use_truth_information ?
            -std::tan( cluster._truth_alpha ):
            -std::tan( cluster._trk_alpha );
          if( std::isnan(talpha) ) continue;

          // z residual
          const auto dz = use_truth_information ? (cluster._z - cluster._truth_z):(cluster._z - cluster._trk_z );
          if( std::isnan(dz) ) continue;
          
          // z angle
          const auto tbeta = use_truth_information ?
            -std::tan( cluster._truth_beta ):
            -std::tan( cluster._trk_beta );
          if( std::isnan(tbeta) ) continue;

          // check against max residuals
          if( std::abs( drp ) > max_residual_drphi ) continue;
          if( std::abs( dz ) > max_residual_dz ) continue;

          // find relevant cell
          auto icell = get_cell( cluster );
          if( icell < 0 ) continue;

          // check pulls
          if( iiter > 0 && nsigmacut[iiter] > 0 )
          {
            
            // find corrections
            const auto drp_corr = get_correction( hphi_rec, cluster );
            const auto dr_corr = get_correction( hr_rec, cluster );
            const auto dz_corr = get_correction( hz_rec, cluster );
            
            const auto pull_rphi = (drp - (drp_corr + talpha*dr_corr))/std::sqrt(erp);
            pulls_drphi[iiter]->Fill( pull_rphi );
            
            const auto pull_z = (dz - (dz_corr + tbeta*dr_corr))/std::sqrt(ez);
            pulls_dz[iiter]->Fill( pull_z );

            if( std::abs(pull_rphi) > nsigmacut[iiter] ) continue;
            if( std::abs(pull_z) > nsigmacut[iiter] ) continue;
          
          }
        
          ++accepted_clusters[iiter];
          ++binentries[icell];
          
          lhs[icell](0,0) += 1./erp;
          lhs[icell](0,1) += 0;
          lhs[icell](0,2) += talpha/erp;
          
          lhs[icell](1,0) += 0;
          lhs[icell](1,1) += 1./ez;
          lhs[icell](1,2) += tbeta/ez;
          
          lhs[icell](2,0) += talpha/erp;
          lhs[icell](2,1) += tbeta/ez;
          lhs[icell](2,2) += square(talpha)/erp + square(tbeta)/ez;
          
          rhs[icell](0,0) += drp/erp;
          rhs[icell](1,0) += dz/ez;
          rhs[icell](2,0) += talpha*drp/erp + tbeta*dz/ez;
          
          // fill histograms
          if( residuals_drphi[icell] ) residuals_drphi[icell]->Fill( drp );
          if( residuals_2d_drphi[icell] ) residuals_2d_drphi[icell]->Fill( talpha, drp );
          if( residuals_dz[icell] ) residuals_dz[icell]->Fill( dz );
          if( residuals_2d_dz[icell] ) residuals_2d_dz[icell]->Fill( tbeta, dz );
          
        }
      }
    }

    // do the inversion
    for( int iphi = 0; iphi < m_phibins; ++iphi )
      for( int ir = 0; ir < m_rbins; ++ir )
      for( int iz = 0; iz < m_zbins; ++iz )
    {
      
      const auto icell = get_cell( iphi, ir, iz );
      
      // cut on number of entries
      // if( binentries[icell]<250 ) continue;
      if( binentries[icell]<10 ) continue;
      
      // save residuals histogram
      if( save_histograms )
      {
        rootFile.Add( residuals_drphi[icell].release() );
        rootFile.Add( residuals_dz[icell].release() );
        if( residuals_2d_drphi[icell] ) rootFile.Add( residuals_2d_drphi[icell].release() );
        if( residuals_2d_dz[icell] ) rootFile.Add( residuals_2d_dz[icell].release() );
      }
      
      // calculate result using linear solving
      const auto cov = lhs[icell].inverse();
      const auto result = lhs[icell].partialPivLu().solve( rhs[icell] );
      
      // fill histograms
      hentries_rec->SetBinContent( iphi+1, ir+1, iz+1, binentries[icell] );
      
      hphi_rec->SetBinContent( iphi+1, ir+1, iz+1, result(0) );
      hphi_rec->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(0,0) ) );
      
      hz_rec->SetBinContent( iphi+1, ir+1, iz+1, result(1) );
      hz_rec->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(1,1) ) );
      
      hr_rec->SetBinContent( iphi+1, ir+1, iz+1, result(2) );
      hr_rec->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(2,2) ) );
      
      out << "Distortions - inverting bin " << iz << ", " << ir << ", " << iphi << std::endl;
      out << "Distortions - entries: " << binentries[icell] << std::endl;
      out << "Distortions - lhs: \n" << lhs[icell] << std::endl;
      out << "Distortions - rhs: \n" << rhs[icell] << std::endl;
      out << "Distortions - drphi: " << result(0) << " +/- " << std::sqrt( cov(0,0) ) << std::endl;
      out << "Distortions - dz: " << result(1) << " +/- " << std::sqrt( cov(1,1) ) << std::endl;
      out << "Distortions - dr: " << result(2) << " +/- " << std::sqrt( cov(2,2) ) << std::endl;
      out << std::endl;
    }

    // print iteration statistics
    std::cout << "Distortions_iterate - iteration: " << iiter << std::endl;
    std::cout << "Distortions_iterate - entries: " << entries << std::endl;
    std::cout << "Distortions_iterate -"
      << " track statistics total: " << total_tracks 
      << " accepted: " << accepted_tracks
      << " fraction: " << 100.*accepted_tracks/total_tracks << "%" 
      << std::endl;
    
    std::cout << "Distortions_iterate -"
      << " cluster statistics total: " << total_clusters 
      << " accepted: " << accepted_clusters[iiter] 
      << " fraction: " << 100.*accepted_clusters[iiter]/total_clusters << "%" 
      << std::endl;

  }
  
  // add non-empty residual histograms  
  for( auto&&h:residuals_drphi ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }
  for( auto&&h:residuals_dz ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }
  for( auto&&h:residuals_2d_drphi ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }
  for( auto&&h:residuals_2d_dz ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }

  // add pulls histograms
  for( auto&&h:pulls_drphi ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }
  for( auto&&h:pulls_dz ) { if( h && h->GetEntries() > 0 ) rootFile.Add( h.release() ); }
  
  // close
  out.close();
  rootFile.Close();

  // print statistics
  std::cout << "Distortions_iterate - entries: " << entries << std::endl;
  std::cout << "Distortions_iterate -"
    << " track statistics total: " << total_tracks 
    << " accepted: " << accepted_tracks
    << " fraction: " << 100.*accepted_tracks/total_tracks << "%" 
    << std::endl;
  
  for( int iiter = 0; iiter < niterations; ++iiter )
  {
    std::cout << "Distortions_iterate -"
      << " iteration: " << iiter 
      << " cluster statistics total: " << total_clusters
      << " accepted: " << accepted_clusters[iiter]
      << " fraction: " << 100.*accepted_clusters[iiter]/total_clusters << "%" 
      << std::endl;
  }
  
  return rootfilename;
}
