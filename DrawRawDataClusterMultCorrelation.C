#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

namespace
{
 
  // detector names
  std::vector<std::string> detector_names;

  MicromegasMapping mapping;
 
  // save detector names
  void save_detector_names()
  {
    
    // get detector names that match tile/layer
    for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name.substr(0,3) ) );
    }
  }

}

//_________________________________________________
TH2* get_scaled_histogram( TH2* source )
{
  const auto name = Form( "%s_scaled", source->GetName() );
  const auto title = Form( "%s (scaled)", source->GetTitle() );
  
  auto dest = new TH2I( name, title, source->GetNbinsX(), 0, 1., source->GetNbinsY(), 0, 1. );
  dest->GetXaxis()->SetTitle( Form( "%s (A.U.)", source->GetXaxis()->GetTitle() ) );
  dest->GetYaxis()->SetTitle( Form( "%s (A.U.)", source->GetYaxis()->GetTitle() ) );
  for( int ix = 0; ix < source->GetNbinsX(); ++ix )
    for( int iy = 0; iy < source->GetNbinsY(); ++iy )
  { dest->SetBinContent( ix+1, iy+1, source->GetBinContent( ix+1, iy+1 ) ); }
  
  dest->GetYaxis()->SetTitleOffset(1.4);
  return dest;
  
}
  
//_________________________________________________
void DrawRawDataClusterMultCorrelation()
{
  set_style(false);
  int runNumber = 9443;  

  const TString inputfilename = Form( "Rootfiles/RawDataClusters-%08i-0000.root", runNumber );
  const TString pdffilename = Form( "Figures/RawDataClusterMultCorrelation-%08i-0000.pdf", runNumber );

  // pdffile
  PdfDocument pdfDocument( pdffilename );

  FileManager fileManager( inputfilename );
  auto h3d = static_cast<TH3*>( fileManager.GetHistogram( "h_cluster_mult_correlation" ) );
  
  save_detector_names();
  
  bool scale_histograms = false;
  
//   if( false )
//   {
//     for( int itile = 0; itile <8; ++itile )
//     {
//       const int detid = itile;
//       auto cvname = Form( "cv_%i", detid);
//       auto cv = std::make_unique<TCanvas>( cvname, cvname, 900,500 );
//       
//       h3d->GetXaxis()->SetRange(detid+1, detid+1);
//       h3d->GetYaxis()->SetRangeUser(1, 50 );
//       h3d->GetZaxis()->SetRangeUser(1, 50 );
// 
//       auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
//       {
//         const auto hname = Form( "h_%s", detector_names[detid].c_str() );
//         const auto htitle = Form( "cluster multiplicity %s (%i)", detector_names[detid].c_str(), runNumber );
//         h2d->SetName( hname );
//         h2d->SetTitle( htitle );
//       }
//         
//       h2d->Draw( "colz" );
//       gPad->SetRightMargin( 0.2 );
//       
//       auto line = new TLine( 0,0,50, 50);
//       line->SetLineStyle( 2 );
//       line->Draw();
//       
//       pdfDocument.Add( cv.get() );
//     }
//   }
  
  if( true )
  {
    auto cv = new TCanvas( "cv", "cv", 1200, 600 );
    cv->Divide( 4, 2 );
    for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile;     
      unsigned int min_cluster = 1;
      unsigned int max_cluster = 30;
      h3d->GetXaxis()->SetRange(detid+1, detid+1);
      h3d->GetYaxis()->SetRangeUser(min_cluster, max_cluster);
      h3d->GetZaxis()->SetRangeUser(min_cluster, max_cluster);

      auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
      {
        const auto hname = Form( "h_%s", detector_names[detid].c_str() );
        // const auto htitle = Form( "Cluster Mult. %s (%i)", detector_names[detid].c_str(), runNumber );
        h2d->SetName( hname );
        h2d->SetTitle( "" );
        h2d->GetXaxis()->SetTitle( "cluster multiplicity #phi" );
        h2d->GetYaxis()->SetTitle( "cluster multiplicity Z" );
        h2d->GetYaxis()->SetTitleOffset(1.3);
      }
      
      cv->cd( detid+1 );
      if( scale_histograms )
      {
        
        get_scaled_histogram(h2d)->Draw("colz");
        auto line = new TLine(0,0,1,1);
        line->SetLineStyle( 2 );
        line->Draw();
        
      } else {
        
        h2d->Draw( "colz" );
        auto line = new TLine(min_cluster, min_cluster, max_cluster, max_cluster);
        line->SetLineStyle( 2 );
        line->Draw();
        
      }

      auto text = new TPaveText(.16,.76,.50,.93, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(1010);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "" );
      text->AddText( "#it{#bf{sPHENIX}} Internal");
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "%s 07/05/2023", detector_names[detid].c_str() ));
      text->Draw();
      
      gPad->SetTopMargin( 0.03 );
      gPad->SetRightMargin( 0.14 );
      gPad->SetLeftMargin( 0.13 );
      
      gPad->SetLogz(true);
      
    }
    pdfDocument.Add( cv );
  }   
   
  if( true )
  {
    auto cv = new TCanvas( "cv1", "cv1", 700, 600 );
    
    unsigned int min_cluster = 1;
    unsigned int max_cluster = 40;
    h3d->GetXaxis()->SetRange(2, 16);
    h3d->GetYaxis()->SetRangeUser(min_cluster, max_cluster);
    h3d->GetZaxis()->SetRangeUser(min_cluster, max_cluster);
    
    auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
    {
      const auto hname = "h_all";
      const auto htitle = Form( "Cluster Mult. (%i)", runNumber );
      h2d->SetName( hname );
      h2d->SetTitle( "" );
      h2d->GetXaxis()->SetTitle( "cluster multiplicity #phi" );
      h2d->GetYaxis()->SetTitle( "cluster multiplicity Z" );
      h2d->GetYaxis()->SetTitleOffset(1.2);
    }
    
    if( scale_histograms )
    {
      
      get_scaled_histogram(h2d)->Draw("colz");
      auto line = new TLine(0,0,1,1);
      line->SetLineStyle( 2 );
      line->Draw();

    } else {
    
      h2d->Draw( "colz" );
      auto line = new TLine(min_cluster, min_cluster, max_cluster, max_cluster);
      line->SetLineStyle( 2 );
      line->Draw();
    
    }
    
    gPad->SetTopMargin( 0.03 );
    gPad->SetRightMargin( 0.16 );
    gPad->SetLeftMargin( 0.14 );

    auto text = new TPaveText(.16,.76,.50,.93, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(1010);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( "" );
    text->AddText( "#it{#bf{sPHENIX}} Internal");
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText( "07/05/2023" );
    text->Draw();
    
    gPad->SetLogz(true);
    pdfDocument.Add( cv );
    
    cv->SaveAs("Figures/TPOT_cluster_multiplicity_correlation.png");
    
  }

}
