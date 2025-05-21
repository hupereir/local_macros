#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TH1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
double get_atan( double y, double x )
{
  const double out = std::atan2( y, x );
  return out > 0 ? out : (out + 2.*M_PI);
}

//____________________________________________________________________________
void Kinematics()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  set_style( false );

  // pdf output
  // const TString tag = "_flat_truth_micromegas_nominal";
  const TString tag = "_flat_truth_micromegas_corrected_mm-coarse_extrapolated-new1";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval_*.root", tag.Data() );
  // const TString inputFile = Form( "DST/CONDOR_%s/dst_eval*.root", tag.Data() );
  const TString pdfFile = Form( "Figures/Kinematics%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  {
    // create canvas
    const TString cvName = "cv_pt";
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    cv->SetRightMargin( 0.2 );

    const TString var = "_tracks._eta:get_atan(_tracks._py,_tracks._px)";
    const TCut cut = "_tracks._ndf > 20";
    auto h = new TH2F( "kinematics", "", 100, 0, 2.*M_PI, 100, -1.5, 1.5 );
    Utils::TreeToHisto( tree, h->GetName(), var, cut, false );

    h->SetTitle( "" );
    h->Draw( "colz" );
    h->GetXaxis()->SetTitle( "#phi (rad)" );
    h->GetYaxis()->SetTitle( "#eta" );

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const double phi = (i+1)*M_PI/6;
      Draw::VerticalLine( cv, phi )->Draw();
    }

    pdfDocument.Add( cv );

    std::cout << "Kinematics -"
      << " tag: " << tag
      << " Entries: " << tree->GetEntries()
      << " track count: " << h->GetEntries()
      << std::endl;

  }

}
