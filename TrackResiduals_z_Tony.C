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

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
float square(float x) { return x*x; }

//____________________________________________________________________________
float get_r( float x, float y ) { return std::sqrt(square(x)+square(y)); }

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}


//_____________________________________________
TString TrackResiduals_z_Tony()
{

  // const TString tag = "_corrected-new";
  const TString tag = "_corrected_notpc-new";
  const TString pdfFile = Form( "Figures/TrackResiduals_z_Tony%s.pdf", tag.Data() );

  std::cout << "TrackResiduals_Tony - pdfFile: " << pdfFile << std::endl;
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager;
  const std::vector<int> runlist = { 41989, 41990, 41991, 41992 };
  for( const auto& runnumber:runlist )
  {
    const TString inputFile = Form( "DST/CONDOR_CombinedDataReconstruction%s/TrackResiduals-%08i-*-full.root", tag.Data(), runnumber );
    std::cout << "TrackResiduals_Tony - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile );
  }

  auto tree = fileManager.GetChain( "residualtree" );

  const TCut trcut(
    "m_crossing == 1"
    "&& m_pt>0.2 && quality<100"
    "&& m_ntpc>20 && m_nmaps>2 && m_nintt>1"
    "&& nmms>0");

  const TCut state_cut("fabs(stategy)<100 && fabs(stategx)<100");
  const TCut layer_cut( "cluslayer>=7&&cluslayer<55" );
  // const TCut layer_cut( "cluslayer>=39&&cluslayer<55" );

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
      const auto hname = Form("htrack_%i", icut );
      auto h = Utils::TreeToHisto( tree, hname, var, trcut&&cuts[icut], true );
      const auto track_count = h->GetEntries();
      std::cout << "TrackResiduals_Tony - icut: " << icut << " cut: " << cuts[icut] << " tracks: " << track_count << std::endl;
    }
  }

  if( true )
  {
    const auto cvname( "cv0" );
    auto cv = new TCanvas( cvname, "", 1400, 1400 );
    cv->Divide( 2, 2 );
    {
      // rphi residuals vs z in TPC
      const TString var = "fabs(clusgr)*(atan2(stategy,stategx)-atan2(clusgy,clusgx))";
      const TString var_2d = Form( "%s:stategz", var.Data() );
      const auto cvname( "cv1" );
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( icut+1 );
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("h_2d_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, -max_dphi, max_dphi );
          Utils::TreeToHisto( tree, hname, var_2d, trcut&&state_cut&&layer_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z (cm)" );
          h->GetYaxis()->SetTitle( "r_{clus}.#Delta#phi_{track-cluster} (cm)" );
          h->Draw( "colz" );
        }

        Draw::PutText( 0.15, 0.85, Form( "Runs 41989-41992 - %s", labels[icut].Data() ) );
        gPad->Update();
        Draw::HorizontalLine(gPad, 0)->Draw();
      }
    }

    {
      // z residuals vs layer
      const TString var = "stategz-clusgz";
      const TString var_2d = Form( "%s:stategz", var.Data() );
      for( int icut = 0; icut < ncuts; ++icut )
      {
        cv->cd( icut+1 + 2 );
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.2);

        {
          const auto hname = Form("h_2d_%i", icut );
          auto h = new TH2F( hname, "", 100, -110, 110, 100, -max_dz, max_dz );
          Utils::TreeToHisto( tree, hname, var_2d, trcut&&layer_cut&&state_cut&&cuts[icut], false );

          h->GetXaxis()->SetTitle( "z (cm)" );
          h->GetYaxis()->SetTitle( "#Deltaz_{track-cluster} (cm)" );
          h->Draw( "colz" );
        }

        Draw::PutText( 0.15, 0.85, Form( "Runs 41989-41992 - %s", labels[icut].Data() ) );
        gPad->Update();
        Draw::HorizontalLine(gPad, 0)->Draw();
      }
    }
    pdfDocument.Add(cv);

  }

  return pdfFile;

}
