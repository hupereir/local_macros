#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>


#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libg4eval_hp.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasDefs.h>

#include <g4eval_hp/MicromegasGeometryContainer.h>

// namespace
// {
//   // geometry container
//   std::unique_ptr<MicromegasGeometryContainer> m_geometry_container;
//
//   // load geometry from file
//   void load_geometry( const std::string& filename )
//   {
//     auto inputfile( TFile::Open( filename.c_str() ) );
//     m_geometry_container.reset( dynamic_cast<MicromegasGeometryContainer*>( inputfile->Get( "MicromegasGeometryContainer" ) ) );
//     m_geometry_container->identify();
//   }
// }

//____________________________________________________________________________
TString make_run_label( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "run %i", runlist[0] );
  return Form( "runs %i-%i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
TString make_run_postfix( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "_%08i", runlist[0] );
  return Form( "_%08i-%08i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//_______________________________________________________________
float square( float value ) { return value*value; }

//_______________________________________________________________
float get_r( float x, float y ) { return std::sqrt( square(x)+square(y) ); }

//_______________________________________________________________
void TpotMatchingEfficiency()
{

  set_style( false );

  MicromegasMapping mapping;


//   // load geometry
//   load_geometry( "micromegas_geometry.root" );

  const TString tag = "_CombinedDataReconstruction_corrected_notpot_masked";
  TString postfix;
  TString run_label;

  FileManager fileManager;
  if( true )
  {
    const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };
    // const std::vector<int> runlist = { 51106, 51107, 51109};
    // const std::vector<int> runlist = { 51249};
    // const std::vector<int> runlist = { 51488};
    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TpcMatchingProfile - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data());
    std::cout << "TpcMatchingProfile - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TpotMatchingEfficiency%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  // file manager
  auto tree = fileManager.GetChain( "T" );

  static constexpr int nbinsx = 20;
  static constexpr int nbinsy = 20;

  static constexpr float max_dx = 20;
  static constexpr float max_dy = 20;

  const TCut trk_cut(
    "_tracks._pt>0.2"
    "&& (_tracks._crossing==1 || _tracks._crossing==0)"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
//     "&& _tracks._nclusters_mvtx>2"
//     "&& _tracks._nclusters_intt>1"
    );


  const double x_min = -10;
  const double x_max = 10;
  const double y_min = -20;
  const double y_max = 20;

  {
    // local coordinates
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->Divide(4, 4);

    // phi layer
    {
      const TString var2d = "_tracks._trk_state_phi._y_local:_tracks._trk_state_phi._x_local";
      const TString var3d = Form( "%s:_tracks._trk_state_phi._tile",var2d.Data() );
      const auto ref_cut = TCut( "_tracks._trk_state_phi._layer>0" ) && trk_cut;
      const TCut max_distance_cut = Form( "fabs(_tracks._best_cluster_phi._x_local-_tracks._trk_state_phi._x_local)<%f", max_dx );
      const TCut found_cut = TCut( "_tracks._best_cluster_phi._layer > 0" ) && max_distance_cut;

      const auto hname_ref( "h_ref_phi");
      auto href_3d = new TH3F( hname_ref, "", 8, 0, 8, nbinsx, -30, 30, nbinsy, -30, 30 );
      Utils::TreeToHisto(tree, href_3d->GetName(), var3d, ref_cut, false);

      const auto hname_found( "h_found_phi" );
      auto hfound_3d = static_cast<TH3*>(href_3d->Clone( hname_found ));
      Utils::TreeToHisto(tree, hfound_3d->GetName(), var3d, ref_cut&&found_cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {

        std::cout << "tile: " << itile << std::endl;
//
//         double x_min = m_geometry_container->get_strip_local_begin(55, itile, 0 ).x();
//         double y_min = m_geometry_container->get_strip_local_begin(55, itile, 0 ).y();
//
//         double x_max = m_geometry_container->get_strip_local_end(55, itile, 255 ).x();
//         double y_max = m_geometry_container->get_strip_local_end(55, itile, 255 ).y();
//
//         if( x_min > x_max ) std::swap( x_min, x_max );
//         if( y_min > y_max ) std::swap( y_min, y_max );
//
//         std::cout << "TpotMatchingEfficiency -"
//           << " x_min: " << x_min
//           << " x_max: " << x_max
//           << " y_min: " << y_min
//           << " y_max: " << y_max
//           << std::endl;

        href_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto href = static_cast<TH2*>(href_3d->Project3D( "zy" ));

        hfound_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto hfound = static_cast<TH2*>(hfound_3d->Project3D( "zy" ));

        const auto hname_eff =  Form( "h_eff_phi_%i", itile);
        auto heff = static_cast<TH2*>(href->Clone( hname_eff ));
        Utils::DivideHistograms2D( hfound, href, heff );
        heff->GetXaxis()->SetTitle( "x_{loc} (cm)" );
        heff->GetYaxis()->SetTitle( "y_{loc} (cm)" );
        heff->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1);
        heff->SetStats(0);
        heff->Draw("colz");
        gPad->SetRightMargin(0.2);
        Draw::PutText( 0.15, 0.85, detname.c_str() );

        auto box = new TBox( x_min, y_min, x_max, y_max );
        box->SetFillStyle(0);
        box->SetLineColor(2);
        box->SetLineStyle(2);
        box->Draw();

        const int ixmin = href->GetXaxis()->FindBin( x_min );
        const int ixmax = href->GetXaxis()->FindBin( x_max );
        const int iymin = href->GetXaxis()->FindBin( y_min );
        const int iymax = href->GetXaxis()->FindBin( y_max );

        const double ref = href->Integral(ixmin, ixmax, iymin, iymax );
        const double found = hfound->Integral(ixmin, ixmax, iymin, iymax );
        const double eff= found/ref;
        std::cout << "TpotMatchingEfficiency - " << detname << " ref: " << ref << " found: " << found << " eff: " << eff << std::endl;
      }

    }

    // z layer
    {
      const TString var2d = "_tracks._trk_state_z._y_local:_tracks._trk_state_z._x_local";
      const TString var3d = Form( "%s:_tracks._trk_state_z._tile",var2d.Data() );
      const auto ref_cut = TCut( "_tracks._trk_state_z._layer>0" ) && trk_cut;
      const TCut max_distance_cut = Form( "fabs(_tracks._best_cluster_z._y_local-_tracks._trk_state_z._y_local)<%f", max_dy );
      const TCut found_cut = TCut( "_tracks._best_cluster_z._layer > 0" ) && max_distance_cut;

      const auto hname_ref( "h_ref_z");
      auto href_3d = new TH3F( hname_ref, "", 8, 0, 8, nbinsx, -30, 30, nbinsy, -30, 30 );
      Utils::TreeToHisto(tree, href_3d->GetName(), var3d, ref_cut, false);

      const auto hname_found( "h_found_z" );
      auto hfound_3d = static_cast<TH3*>(href_3d->Clone( hname_found ));
      Utils::TreeToHisto(tree, hfound_3d->GetName(), var3d, ref_cut&&found_cut, false);

      // loop over tiles
      for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
      {
//
//         double x_min = m_geometry_container->get_strip_local_begin(55, itile, 0 ).x();
//         double y_min = m_geometry_container->get_strip_local_begin(55, itile, 0 ).y();
//
//         double x_max = m_geometry_container->get_strip_local_end(55, itile, 255 ).x();
//         double y_max = m_geometry_container->get_strip_local_end(55, itile, 255 ).y();
//         if( x_min > x_max ) std::swap( x_min, x_max );
//         if( y_min > y_max ) std::swap( y_min, y_max );

        href_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto href = static_cast<TH2*>(href_3d->Project3D( "zy" ));

        hfound_3d->GetXaxis()->SetRange( itile+1, itile+1 );
        auto hfound = static_cast<TH2*>(hfound_3d->Project3D( "zy" ));

        const auto hname_eff =  Form( "h_eff_z_%i", itile);
        auto heff = static_cast<TH2*>(href->Clone( hname_eff ));
        Utils::DivideHistograms2D( hfound, href, heff );
        heff->GetXaxis()->SetTitle( "x_{loc} (cm)" );
        heff->GetYaxis()->SetTitle( "y_{loc} (cm)" );
        heff->SetTitle("");

        // get detector name
        const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
        const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

        cv->cd(2*itile+1+1);
        heff->SetStats(0);
        heff->Draw("colz");
        gPad->SetRightMargin(0.2);
        Draw::PutText( 0.15, 0.85, detname.c_str() );

        auto box = new TBox( x_min, y_min, x_max, y_max );
        box->SetFillStyle(0);
        box->SetLineColor(2);
        box->SetLineStyle(2);
        box->Draw();

        const int ixmin = href->GetXaxis()->FindBin( x_min );
        const int ixmax = href->GetXaxis()->FindBin( x_max );
        const int iymin = href->GetXaxis()->FindBin( y_min );
        const int iymax = href->GetXaxis()->FindBin( y_max );

        const double ref = href->Integral(ixmin, ixmax, iymin, iymax );
        const double found = hfound->Integral(ixmin, ixmax, iymin, iymax );
        const double eff= found/ref;
        std::cout << "TpotMatchingEfficiency - " << detname << " ref: " << ref << " found: " << found << " eff: " << eff << std::endl;

      }

    }

    pdfDocument.Add(cv);
  }


}
