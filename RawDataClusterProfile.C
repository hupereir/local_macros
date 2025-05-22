
#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

//_____________________________________________________________________________
void RawDataClusterProfile(int runNumber = 21082)
// void RawDataClusterProfile(int runNumber = 9467)
{

  set_style( false );

  MicromegasMapping mapping;

  const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/dst_eval_masked-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataClusterProfile-%08i-0000.pdf", runNumber );

  std::cout << "RawDataClusterProfile - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataClusterProfile - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  const TString var2dx = "clusters[].x_local:(clusters[].tile+8*(clusters[].layer-55))";
  const TString var2dy = "clusters[].y_local:(clusters[].tile+8*(clusters[].layer-55))";
  const auto h2dx = new TH2I( "h2dx", "", 16, 0, 16, 100, -15, 15 );
  const auto h2dy = new TH2I( "h2dy", "", 16, 0, 16, 100, -30, 30 );
  for( const auto& h: {h2dx, h2dy} )
  {
    h->StatOverflows(true);
    h->GetYaxis()->SetTitle( "x_{local} (cm)" );
    h->GetZaxis()->SetTitle( "A.U." );
  }

  const TCut full_event_cut( "_events[0]._nrawhits_micromegas>4000" );
  Utils::TreeToHisto( tree, h2dx->GetName(), var2dx, full_event_cut, false );
  Utils::TreeToHisto( tree, h2dy->GetName(), var2dy, full_event_cut, false );

  const auto h_phi = new TH1F( "h_phi", "", 100, -15, 15 );
  const auto h_z = new TH1F( "h_z", "", 100, -30, 30 );
  for( auto&& h:{h_phi,h_z})
  {
    h->StatOverflows(true);
    h->GetXaxis()->SetTitle( "x_{local} (cm)" );
    h->GetYaxis()->SetTitle( "A.U." );
    h->GetYaxis()->SetTitleOffset( 1.8 );
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );
  }

  const std::vector<std::string> detectors_phi =
  {
    "SCIP",
    "NCOP",
    "SEIP",
    "NEIP",
    "SWIP",
    "NWIP"
  };


  const std::vector<std::string> detectors_z =
  {
    "SCIZ",
    "NCIZ",
    "NCOZ",
    "SEIZ",
    "SWIZ",
    "NWIZ"
  };

  // create canvas and divide
  const auto cvname = "cv";
  auto cv = new TCanvas( cvname, cvname, 900, 900 );
  cv->Divide(4,4);

  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer <2; ++ilayer )
  {

    // get actual layer and segmentation
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    const bool is_phi = (layer == 55);

    for( int itile = 0; itile <8; ++itile )
    {

      // generate hitset key from layer, tile and segmentation. This allows to retrieve the detector name from the Mapping
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );

      const auto hname = Form( "h_%i_%i", ilayer, itile );
      const auto htitle = Form( "%s - %i,%i", name.c_str(), layer, itile );

      // get the bin matching layer, tile and region in the 3D histgoram
      const int bin = itile + 8*ilayer;
      const auto h2d = (layer == 55) ? h2dx:h2dy;
      h2d->GetXaxis()->SetRange( bin+1, bin+1 );

      auto h = static_cast<TH1*>( h2d->ProjectionY() );
      h->SetName( hname );
      h->SetTitle( htitle );

      // sum to average
      if( is_phi )
      {
        if( std::find( detectors_phi.begin(), detectors_phi.end(), name ) != detectors_phi.end() )
        { h_phi->Add( h ); }

      } else {

        if( std::find( detectors_z.begin(), detectors_z.end(), name ) != detectors_z.end() )
        { h_z->Add( h ); }

      }

      // draw
      cv->cd(bin+1);
      h->Scale( 1./h->GetEntries() );
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow );
      h->Draw( "hist" );

      // get the mean
      const double mean = double(h->GetEntries())/h->GetXaxis()->GetNbins();

      {
        auto text = new TPaveText(0.37,0.76,0.94,0.9, "NDC" );
        text->SetFillColor(0);
        text->SetFillStyle(0);
        text->SetBorderSize(0);
        text->SetTextAlign(11);
        text->AddText(Form( "Mean: %.2f", mean ));
        text->Draw();
      }

    }

  }

  pdfDocument.Add( cv );


  {
    auto cv = new TCanvas( "cvp", "cv", 980, 900 );
    h_phi->Scale( 1./h_phi->GetEntries() );
    h_phi->Draw("hist");
    gPad->SetLeftMargin( 0.18 );

    const double n_cluster_mean = h_phi->GetMean();
    {
      auto text = new TPaveText(0.3,0.18,0.92,0.31, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "#it{#bf{sPHENIX}}" );
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText("#phi strips");
      text->Draw();
    }
    pdfDocument.Add(cv);

    // cv->SaveAs( Form( "Figures/RawDataClusterProfile_phi-%08i-0000.png", runNumber ) );
    cv->SaveAs( Form( "Figures/RawDataClusterProfile_phi-%08i-0000.pdf", runNumber ) );
  }

  {
    auto cv = new TCanvas( "cvz", "cv", 980, 900 );
    h_z->Scale( 1./h_z->GetEntries() );
    h_z->Draw("hist");
    gPad->SetLeftMargin( 0.18 );

    {
      auto text = new TPaveText(0.3,0.18,0.92,0.31, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "#it{#bf{sPHENIX}}" );
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText("z strips");
      text->Draw();
    }
    pdfDocument.Add(cv);
    // cv->SaveAs( Form( "Figures/RawDataClusterProfile_z-%08i-0000.png", runNumber ) );
    cv->SaveAs( Form( "Figures/RawDataClusterProfile_z-%08i-0000.pdf", runNumber ) );
  }

  return pdfFile;

}
