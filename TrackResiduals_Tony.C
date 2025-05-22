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
TString TrackResiduals_Tony()
{

  // const TString tag = "_flat_acts_truth_notpc_nodistortion";
  // const TString tag = "_corrected-new";
  const TString tag = "_uncorrected-0";
  // const TString tag = "_corrected-1";
  TString postfix;
  TString run_label;

  FileManager fileManager;

  if( true )
  {
    // const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };
    // const std::vector<int> runlist = {49714, 49715, 49716, 49717, 49718, 49719, 49720, 49721, 49722 };
    // const std::vector<int> runlist = {50015, 50016, 50017, 50018, 50019, 50020};
    // const std::vector<int> runlist = { 50450, 50451, 50452, 50454 };
    // const std::vector<int> runlist = { 51249 };
    // const std::vector<int> runlist = { 51488 };
    // const std::vector<int> runlist = { 53199 };
    const std::vector<int> runlist = { 53285 };

    for( const auto& runnumber:runlist )
    {
      const TString inputFile = Form( "DST/CONDOR_CombinedDataReconstruction%s/TrackResiduals-%08i-*-full.root", tag.Data(), runnumber );
      std::cout << "TrackResiduals_Tony - inputFile: " << inputFile << std::endl;
      fileManager.AddFiles( inputFile );
    }

    run_label = make_run_label( runlist );
    postfix = make_run_postfix( runlist );

  } else {

    const TString inputFile = Form( "DST/CONDOR%s/TrackResiduals_flat_acts_truth_notpc_nodistortion_?.root", tag.Data() );
    std::cout << "TrackResiduals_Tony - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );

  }

  const TString pdfFile = Form( "Figures/TrackResiduals_Tony%s%s.pdf", tag.Data(),postfix.Data());
  const TString rootFilename = Form( "Rootfiles/TrackResiduals_Tony%s%s.root", tag.Data(),postfix.Data());

  std::cout << "TrackResiduals_Tony - pdfFile: " << pdfFile << std::endl;
  std::cout << "TrackResiduals_Tony - rootFilename: " << rootFilename << std::endl;
  RootFile rootFile( rootFilename );
  PdfDocument pdfDocument( pdfFile );

  auto tree = fileManager.GetChain( "residualtree" );
  std::cout << "entries: " << tree->GetEntries() << std::endl;
  const TCut trcut(
    "(m_crossing == 1 || m_crossing == 0)"
    "&& m_pt>0.2 "
    "&& quality<100"
    "&& m_ntpc>20"
    "&& m_nmaps>2"
    "&& m_nintt>1"
    // "&& nmms>0"
    );

  const TCut state_cut("fabs(stategy)<100 && fabs(stategx)<100");

  const TCut cuts[] =
  {
    "charge<0",
    "charge>0"
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
      const TString var = "crossing";
      const auto hname = Form("hcrossing_%i", icut );
      auto h = Utils::TreeToHisto( tree, hname, var, trcut&&cuts[icut], true );
      const auto track_count = h->GetEntries();
      std::cout << "TrackResiduals_Tony - icut: " << icut << " cut: " << cuts[icut] << " tracks: " << track_count << std::endl;

      rootFile.Add( h );
    }
  }

  if( true )
  {
    auto cv = new TCanvas( "cv0", "", 1400, 700 );
    cv->Divide( 2, 1 );

    // track pseudo rapidity measurement
    const TString var = "get_eta(atan2(m_pt,m_pz))";

    for( int icut = 0; icut < ncuts; ++icut )
    {
      cv->cd( icut+1 );
      gPad->SetTopMargin(0.05);

      const auto hname = Form("heta_%i", icut );
      auto h = new TH1F( hname, "", 100, -1.5, 1.5 );
      Utils::TreeToHisto( tree, hname, var, trcut&&cuts[icut], false );
      h->GetXaxis()->SetTitle("#eta");

      h->Draw();
      rootFile.Add( h );
    }
    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto cv = new TCanvas( "cv1", "", 1400, 1400 );
    cv->Divide( 2, 2 );

    {
      // clusters on track
      // x,y plots

      const TString var = "clusgy:clusgx";
      const TString var2 = "stategy:stategx";

      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("hxy_cluster_%i", icut );
          auto h = new TH2F( hname, "", 170, -90, 90, 170, -90, 90 );
          Utils::TreeToHisto( tree, hname, var, trcut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "x_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "y_{clus} (cm)" );
          h->Draw();
          rootFile.Add( h );
        }

        {
          const auto hname = Form("hxy_state_%i", icut );
          auto h = new TH2F( hname, "", 170, -90, 90, 170, -90, 90 );
          Utils::TreeToHisto( tree, hname, var2, trcut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "x_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "y_{clus} (cm)" );
          h->SetLineColor(2);
          h->SetMarkerColor(2);
          h->Draw("same");
          rootFile.Add( h );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }

    // clusters on track
    // r,z plots
    {
      const TString var = "get_r(clusgy,clusgx):clusgz";
      const TString var2 = "get_r(stategy,stategx):stategz";
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( 2+ icut+1 );
        gPad->SetTopMargin(0.05);

        {
          const auto hname = Form("hrz_cluster_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, 0, 90 );
          Utils::TreeToHisto( tree, hname, var, trcut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus} (cm)" );
          h->Draw();
          rootFile.Add( h );
        }

        {
          const auto hname = Form("hrz_state_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, 0, 90 );
          Utils::TreeToHisto( tree, hname, var2, trcut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z_{clus} (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus} (cm)" );
          h->SetLineColor(2);
          h->SetMarkerColor(2);
          h->Draw("same");
          rootFile.Add( h );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
      }
    }
    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto cv = new TCanvas( "cv2", "", 1400, 1400 );
    cv->Divide( 2, 2 );
    {
      // rphi residuals vs layer
      const TString var = "fabs(clusgr)*(atan2(stategy,stategx)-atan2(clusgy,clusgx))";
      const TString var_2d = Form( "%s:cluslayer", var.Data() );
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("hdrphi_layer_%i", icut );
          auto h = new TH2F( hname, "", 60, 0, 60, 100, -max_dphi, max_dphi );
          Utils::TreeToHisto( tree, hname, var_2d, trcut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "layer" );
          h->GetYaxis()->SetTitle( "r_{clus}.#Delta#phi_{track-cluster} (cm)" );
          h->Draw( "colz" );
          rootFile.Add( h );
        }

        Draw::PutText( 0.15, 0.85, Form( "%s- %s", run_label.Data(), labels[icut].Data() ) );
        gPad->Update();
        Draw::HorizontalLine(gPad, 0)->Draw();
      }
    }

    {
      // z residuals vs layer
      const TString var = "stategz-clusgz";
      const TString var_2d = Form( "%s:cluslayer", var.Data() );
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( icut+1 + 2 );
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("hdz_layer_%i", icut );
          auto h = new TH2F( hname, "", 60, 0, 60, 100, -max_dz, max_dz );
          Utils::TreeToHisto( tree, hname, var_2d, trcut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "layer" );
          h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );
          h->Draw( "colz" );
          rootFile.Add( h );
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
