#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <RootUtil/RootFile.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
float square(float x) { return x*x; }

//____________________________________________________________________________
float get_r( float x, float y ) { return std::sqrt(square(x)+square(y)); }

//____________________________________________________________________________
float get_eta( float theta )
{ return -std::log(std::tan(theta/2)); }

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

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

//_____________________________________________
TString TrackResiduals_hp()
{


  // const TString tag = "_flat_acts_truth_notpc_nodistortion";
  // const TString tag = "_uncorrected-2";
  // const TString tag = "_corrected-2";
  // const TString tag = "_corrected_module_edge-2";
  const TString tag = "_CombinedDataReconstruction_zf_corrected";
  TString postfix;
  TString run_label;

  FileManager fileManager;

  {
    // zero field runs
    const std::vector<int> runlist = { 52077, 52078 };
    // const std::vector<int> runlist = { 52401 };
    // const std::vector<int> runlist = { 53199 };
    // const std::vector<int> runlist = { 53285 };

    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR%s/dst_eval-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TrackResiduals_hp - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );
  }

  const bool match_silicons = false;
  const TString pdfFile = match_silicons ?
    Form( "Figures/TrackResiduals_hp_matched%s%s.pdf", tag.Data(),postfix.Data()):
    Form( "Figures/TrackResiduals_hp%s%s.pdf", tag.Data(),postfix.Data());
  const TString rootFilename = match_silicons ?
    Form( "Rootfiles/TrackResiduals_hp_matched%s%s.root", tag.Data(),postfix.Data()):
    Form( "Rootfiles/TrackResiduals_hp%s%s.root", tag.Data(),postfix.Data());

  std::cout << "TrackResiduals_hp - pdfFile: " << pdfFile << std::endl;
  std::cout << "TrackResiduals_hp - rootFilename: " << rootFilename << std::endl;
  RootFile rootFile( rootFilename );
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "T" );

  std::cout << "TrackResiduals_hp - entries: " << tree->GetEntries() << std::endl;


  const TCut base_trk_cut(
    "(_tracks._crossing == 1 || _tracks._crossing == 0)"
    "&& _tracks._pt>0.3"
    "&& _tracks._ndf > 0"
    "&& (_tracks._chisquare/_tracks._ndf)<100"
    "&& _tracks._nclusters_tpc>20"
    );

  const TCut match_silicon_cut(
    "_tracks._nclusters_mvtx>=2"
    "&& _tracks._nclusters_intt>=1"
    );

  const TCut trk_cut = match_silicons ?
    base_trk_cut && match_silicon_cut:
    base_trk_cut;

  const TCut state_cut(
    "_tracks._clusters._trk_phi_error>0"
    "&&  _tracks._clusters._trk_z_error>0");


  const TCut cuts[] =
  {
    "_tracks._charge<0",
    "_tracks._charge>0"
  };

  const TString labels[] =
  {
    "charge<0",
    "charge>0"
  };

  static constexpr int ncuts = 2;

//   const double max_dphi = 0.5;
//   const double max_dz = 2;

  const double max_dphi = 2;
  const double max_dz = 5;

  // track count
  if( true )
  {
    for( int icut = 0; icut < ncuts; ++icut )
    {
      const TString var = "_tracks._crossing";
      const auto hname = Form("htrack_%i", icut );
      auto h = Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], true );
      const auto track_count = h->GetEntries();
      std::cout << "TrackResiduals_hp - icut: " << icut << " cut: " << cuts[icut] << " tracks: " << track_count << std::endl;
    }
  }

  if( true )
  {
    auto cv = new TCanvas( "cv0", "", 1400, 700 );
    cv->Divide( 2, 1 );

    // track pseudo rapidity measurement
    const TString var = "get_eta(atan2(_tracks._pt,_tracks._pz))";

    for( int icut = 0; icut < ncuts; ++icut )
    {
      cv->cd( icut+1 );
      gPad->SetTopMargin(0.05);

      const auto hname = Form("heta_%i", icut );
      auto h = new TH1F( hname, "", 100, -1.5, 1.5 );
      Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], false );
      h->GetXaxis()->SetTitle("#eta");

      h->Draw();
      rootFile.Add( h );
    }
    pdfDocument.Add(cv);
  }

  if( false )
  {
    const auto cvname( "cv0" );
    auto cv = new TCanvas( cvname, "", 1400, 1400 );
    cv->Divide( 2, 2 );

    {
      // clusters on track
      // x,y plots
      const TString var = "_tracks._clusters._y:_tracks._clusters._x";
      const TString var2 = "_tracks._clusters._trk_y:_tracks._clusters._trk_x";

      for( int icut = 0; icut < ncuts; ++icut )
      {

        std::cout << "TrackingResiduals_hp - icut: " << icut << std::endl;

        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("h_%i", icut );
          auto h = new TH2F( hname, "", 170, -85, 85, 170, -85, 85 );
          Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "x_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "y_{clus} (cm)" );
          h->Draw();
        }

        {
          const auto hname = Form("h2_%i", icut );
          auto h = new TH2F( hname, "", 170, -85, 85, 170, -85, 85 );
          Utils::TreeToHisto( tree, hname, var2, trk_cut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "x_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "y_{clus} (cm)" );
          h->SetLineColor(2);
          h->SetMarkerColor(2);
          h->Draw("same");
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }

    // clusters on track
    // r,z plots
    {
      const TString var = "_tracks._clusters._r:_tracks._clusters._z";
      const TString var2 = "_tracks._clusters._trk_r:_tracks._clusters._trk_z";
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( 2+ icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("hrz_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, 0, 90 );
          Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus} (cm)" );
          h->Draw();
        }

        {
          const auto hname = Form("hrz2_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, 0, 90 );
          Utils::TreeToHisto( tree, hname, var2, trk_cut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus} (cm)" );
          h->SetLineColor(2);
          h->SetMarkerColor(2);
          h->Draw("same");
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }
    pdfDocument.Add(cv);
  }

  if( true )
  {
    const auto cvname( "cv0" );
    auto cv = new TCanvas( cvname, "", 1400, 1400 );
    cv->Divide( 2, 2 );

    {
      // clusters on track
      // x,y plots
      const TString var = "_tracks._clusters._y:_tracks._clusters._x";
      for( int icut = 0; icut < ncuts; ++icut )
      {

        std::cout << "TrackingResiduals_hp - icut: " << icut << std::endl;

        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("h_%i", icut );
          auto h = new TH2F( hname, "", 170, -85, 85, 170, -85, 85 );
          Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "x_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "y_{clus} (cm)" );
          h->Draw("col");
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }

    // clusters on track
    // r,z plots
    {
      const TString var = "_tracks._clusters._r:_tracks._clusters._z";
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( 2+ icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("hrz_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, 0, 90 );
          Utils::TreeToHisto( tree, hname, var, trk_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus} (cm)" );
          h->Draw( "col" );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }
    pdfDocument.Add(cv);
  }

  if( true )
  {
    const auto cvname( "cv1" );
    auto cv = new TCanvas( cvname, "", 1400, 1400 );
    cv->Divide( 2, 2 );
    {
      // rphi residuals vs layer
      const TString var = "_tracks._clusters._r*delta_phi(_tracks._clusters._trk_phi-_tracks._clusters._phi)";
      const TString var_2d = Form( "%s:_tracks._clusters._layer", var.Data() );

      for( int icut = 0; icut < ncuts; ++icut )
      {

        std::cout << "TrackingResiduals_hp - icut: " << icut << std::endl;

        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("h_2d_%i", icut );
          auto h = new TH2F( hname, "", 60, 0, 60, 100, -max_dphi, max_dphi );
          Utils::TreeToHisto( tree, hname, var_2d, trk_cut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "layer" );
          h->GetYaxis()->SetTitle( "r_{clus}.#Delta#phi_{track-cluster} (cm)" );
          h->Draw( "colz" );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
        Draw::HorizontalLine(gPad, 0)->Draw();
      }
    }

    {
      // z residuals vs layer
      const TString var = "_tracks._clusters._trk_z-_tracks._clusters._z";
      const TString var_2d = Form( "%s:_tracks._clusters._layer", var.Data() );

      for( int icut = 0; icut < ncuts; ++icut )
      {
        std::cout << "TrackingResiduals_hp - icut: " << icut << std::endl;

        cv->cd( icut+1 + 2);
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("h_2d_%i", icut );
          auto h = new TH2F( hname, "", 60, 0, 60, 100, -max_dz, max_dz);
          Utils::TreeToHisto( tree, hname, var_2d, trk_cut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "layer" );
          h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );
          h->Draw( "colz" );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
        Draw::HorizontalLine(gPad, 0)->Draw();
      }

    }
    pdfDocument.Add(cv);

  }

  return pdfFile;

}
