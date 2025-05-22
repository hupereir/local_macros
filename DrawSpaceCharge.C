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
#include <THnSparse.h>
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//________________________________________________________________
TString DrawSpaceCharge( TString tag = TString() )
{

  set_style( false );

  if( tag.IsNull() ) tag = "_new-ibf03_2";
  const TString inputFile = Form( "Rootfiles/spacechargemap_integrated%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/spacechargemap_integrated%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // open TFile
  TFile f( inputFile );
  if( !f.IsOpen() ) return TString();

  // get histogram and histogram dimentsions
  auto hn = static_cast<THnSparse*>( f.Get( "integral" ) );
  std::array<int,3> bins;
  for( int i = 0; i < 3; i++ )
  {
    auto axis = hn->GetAxis(i);
    bins[i] = axis->GetNbins();
  }

  // reduce axis range to avoid edge effects
  hn->GetAxis(0)->SetRange( 1, bins[0]-1 );
  hn->GetAxis(1)->SetRange( 2, bins[1]-1 );

  {
    // phi vs z
    const TString cvName( "cv_phiz" );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    cv->SetRightMargin( 0.2 );

    auto h = hn->Projection( 2, 0 )->RebinX(10);
    h->SetTitle( "" );
    h->Draw( "colz" );
    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );
  }

  {
    // r vs z
    const TString cvName( "cv_rz" );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    cv->SetRightMargin( 0.2 );
    auto h = hn->Projection( 1, 0 )->RebinX(10);
    h->SetTitle( "" );
    h->Draw( "colz" );
    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );
  }

  {
    // z
    const TString cvName( "cv_z" );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    auto h = hn->Projection( 0 )->RebinX(10);
    h->SetTitle( "" );
    h->SetMinimum(0);
    h->Draw();
    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );
  }

  return pdfFile;

}
