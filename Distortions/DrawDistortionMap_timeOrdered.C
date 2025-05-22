#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

TString DrawDistortionMap_timeOrdered()
{

  set_style( false );

  const TString tag = "TimeOrderedDistortions_deltas-new";
  const TString inputfile = "/sphenix/user/rcorliss/distortion_maps/2022.07/TimeOrderedDistortions_deltas.root";

  // open TFile
  auto f = TFile::Open( inputfile );
  if( !f ) return TString();
  
  // tree  
  auto tree = static_cast<TTree*>( f->Get("TimeDists") );
  if( !tree ) 
  {
    std::cout << "ConvertDistortions_timeOrdered - invalid tree" << std::endl;
    return TString();
  }
  
  TString pdfFile( Form( "Figures/DistortionMap%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  const TString prefix[] = {
    "hIntDistortionR",
    "hIntDistortionP",
    "hIntDistortionZ" 
  };
  
  const TString suffix[] = {
    "_negz",
    "_posz" 
  };

  const TString label[] = {
    "#Deltar (cm)",
    "r#Delta#phi (cm)",
    "#Deltaz (cm)"
  };

  // set branches
  using array_t = std::array<TH3*, 2>;  
  array_t histograms[3];
  for( int isuffix = 0; isuffix < 2; ++isuffix )
    for( int i = 0; i<3; ++i )
  {
    histograms[i][isuffix] = new TH3F;
    const auto branchName = Form( "%s%s", prefix[i].Data(), suffix[isuffix].Data() );
    std::cout << "DrawDistortionMap_timeOrdered - branchName: " << branchName << std::endl;
    tree->SetBranchAddress( branchName, &((histograms[i])[isuffix]) );
  }
  
  // loop over entries
  const auto entries = tree->GetEntries();
  // const int entries = 10;
  std::cout << "DrawDistortionMap_timeOrdered - entries: " << entries << std::endl;
  for( int ientry = 0; ientry < entries; ++ientry )
  {
    if( !tree->GetEntry(ientry) ) continue;
    std::cout << "DrawDistortionMap_timeOrdered - entry: " << ientry << std::endl;
  
    const TString cvname = Form( "cv_%i", ientry );
    TCanvas* cv( new TCanvas( cvname, cvname, 1200, 800 ) );
    cv->Divide( 3, 2 );
  
    for( int isuffix = 0; isuffix < 2; ++isuffix )
    {
      
      for( int ih = 0; ih < 3; ++ih )
      {
        auto h3= histograms[ih][isuffix];
        
        std::cout << "DrawDistortionMap_timeOrdered - histogram: " << h3->GetName() << std::endl;
        
        h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
        auto h = h3->Project3D( "yz" );
        h->SetTitle( "" );
        
        h->Scale( 1./h3->GetXaxis()->GetNbins() );
        h->GetXaxis()->SetTitle( "z (cm)" );
        h->GetYaxis()->SetTitle( "r (cm)" );
        h->GetZaxis()->SetTitle( label[ih] );
        h->GetZaxis()->SetTitleOffset( 1.7 );

        cv->cd(ih+3*isuffix+1);
        gPad->SetRightMargin( 0.22 );
      
        h->Draw("colz");
        Draw::HorizontalLine( gPad, 30 )->Draw();

        std::cout << "maximum: " << h->GetMaximum() << std::endl;
        std::cout << "minimum: " << h->GetMinimum() << std::endl;
        
      }
    }

    pdfDocument.Add( cv );
  }
  
  return pdfFile;

}
