#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <micromegas/MicromegasMapping.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

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

//_____________________________________________________
void RawDataClusterCountCorrelation_new()
{
  const TString inputfilename = "DST/CONDOR_CombinedDataEvaluation/dst_eval-00020121-0000-full.root";
  const TString pdffilename = "Figures/RawDataClusterCountCorrelation_new-00020121-0000-full.pdf";
  
  FileManager f(inputfilename);
  auto tree = f.GetChain( "T" );
  
  PdfDocument pdfDocument( pdffilename );
  
  save_detector_names();

  if( true )
  {
    auto cv = new TCanvas( "cv0", "cv0", 980, 900 );
    const auto var2d = "n_z_clusters:n_phi_clusters";
    const auto hname = "hall";
    auto h = new TH2F( hname, hname, 100, 0, 200, 100, 0, 200 );
    h->GetXaxis()->SetTitle( "N_{clusters} #phi" );
    h->GetYaxis()->SetTitle( "N_{clusters} z" );
    h->GetYaxis()->SetTitleOffset(1.3);

    Utils::TreeToHisto( tree, hname, var2d, TCut(), false );
    h->SetTitle( "" );
   
    h->Draw( "colz" );
 
    gPad->SetTopMargin( 0.03 );
    gPad->SetRightMargin( 0.14 );
    gPad->SetLeftMargin( 0.13 );
    
    auto text = new TPaveText(.17,.86,.55,.94, "NDC" );
    text->SetFillColor(0);
    // text->SetFillStyle(1010);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText( "" );
    // text->AddText( "#it{#bf{sPHENIX}} Internal");
    // text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText( "RUN 20121" );
    text->Draw();
    
    gPad->SetLogz(true);    
    pdfDocument.Add( cv );
  }
  

  if( true )
  {
    auto cv = new TCanvas( "cv", "cv", 1200, 600 );
    cv->Divide( 4, 2 );
    for( int itile = 0; itile <8; ++itile )
    {
      const auto var2d = Form( "n_detector_clusters[%i]:n_detector_clusters[%i]", itile+8, itile );
      const auto hname = Form( "h%i", itile );
      auto h = new TH2F( hname, hname, 30, 0, 30, 30, 0, 30 );
      h->GetXaxis()->SetTitle( "N_{clusters} #phi" );
      h->GetYaxis()->SetTitle( "N_{clusters} z" );
      h->GetYaxis()->SetTitleOffset(1.3);
      
      Utils::TreeToHisto( tree, hname, var2d, TCut(), false );
      h->SetTitle( "" );

      cv->cd( itile+1 );
      h->Draw( "colz" );
 
      gPad->SetTopMargin( 0.03 );
      gPad->SetRightMargin( 0.14 );
      gPad->SetLeftMargin( 0.13 );

      auto text = new TPaveText(.16,.71,.62,.93, "NDC" );
      text->SetFillColor(0);
      // text->SetFillStyle(1010);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "" );
      // text->AddText( "#it{#bf{sPHENIX}} Internal");
      // text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText( "RUN 20121" );
      text->AddText( detector_names[itile].substr(0,3).c_str() );
      text->Draw();
      
      
      gPad->SetLogz(true);    
      
    }
    pdfDocument.Add( cv );
  }  
  
}
