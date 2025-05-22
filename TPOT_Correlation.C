#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

//____________________________________________________________________________
class FullEventContainer: public TObject
{
  
  public:
  
  void Reset()
  {
    lvl1_bco = 0;
    lvl1_counter = 0;
    n_waveforms_all = 0;
    n_waveforms_signal = 0;
    n_clusters = 0;
    n_phi_clusters = 0;
    n_z_clusters = 0;
    n_detector_clusters.clear();
  }
  
  /// ll1 bco
  uint64_t lvl1_bco = 0;

  /// ll1 counter
  uint32_t lvl1_counter = 0;
  
  // number of waveforms
  unsigned short n_waveforms_all = 0;
  
  // number of signal waveforms
  unsigned short n_waveforms_signal = 0;
  
  // number of clusters
  unsigned short n_clusters = 0;
  
  // number of clusters per detector
  std::vector<unsigned short> n_detector_clusters;
  
  // number of clusters
  unsigned short n_phi_clusters = 0;
  
  // number of clusters
  unsigned short n_z_clusters = 0;
    
  ClassDef(FullEventContainer,1)
  
};

namespace
{

  std::vector<std::string> module_names;
  std::vector<std::string> detector_names;
  MicromegasMapping mapping;

  // save detector names
  void save_detector_names()
  {
    
    // get detector names that match tile/layer
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( name );
      
      // module name
      if( ilayer == 0 ) module_names.push_back( name.substr(0,3) );
      
    }
  }
  
}

//_________________________________________________
void TPOT_Correlation( const int runnumber = 20445 )
{
  set_style( false );
  save_detector_names();
  
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  const TString pdffilename = Form( "Figures/TPOT_Correlation-%08i-0000.pdf", runnumber );
  const TString pngfilename = Form( "Figures/TPOT_Correlation-%08i-0000.png", runnumber );
  const TString rootfilename = Form( "Rootfiles/TPOT_Correlation-%08i-0000.root", runnumber );

  std::cout << "TPOT_Correlation - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  PdfDocument pdfDocument( pdffilename );
  RootFile rootfile( rootfilename );
  
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  auto tpot_container = new FullEventContainer;
  tpot_tree->SetBranchAddress( "Event", &tpot_container );

  // correlation histogram
  int ioffset = 0;
  const auto h_name = "h_correlation";
  auto h_correlation = new TH2F( h_name, "", 180, 0, 180, 180, 0, 180 );
  h_correlation->GetXaxis()->SetTitle( "TPOT N_{clusters} #phi" );
  h_correlation->GetYaxis()->SetTitle( "TPOT N_{clusters} z" );
  rootfile.Add( h_correlation );

  std::array<TH1*,8> h_module_correlation;
  for( int tile = 0; tile < 8; ++tile )
  {
    const auto h_name = Form( "h_correlation_%i", tile );
    const auto h_title = Form( "TPOT correlation %s - offset=%i", module_names[tile].c_str(), ioffset );
    h_module_correlation[tile] = new TH2F( h_name, h_title, 30, 0, 30, 30, 0, 30 );
    h_module_correlation[tile]->GetXaxis()->SetTitle( "N_{clusters} #phi" );
    h_module_correlation[tile]->GetYaxis()->SetTitle( "N_{clusters} z" );
    rootfile.Add( h_module_correlation[tile] );
  }
  
  const int tpot_entries = tpot_tree->GetEntries();
  std::cout << "TPOT_Correlation -"
    << " tpot_entries: " << tpot_entries
    << std::endl;
  
  for( int tpot_entry = 0; tpot_entry < tpot_entries; ++tpot_entry )
  {
  
    if( !(tpot_entry%100) ) std::cout << "TPOT_Correlation - entry: " << tpot_entry << std::endl;
    
    // get tpot entry
    tpot_tree->GetEntry(tpot_entry);
    // const auto n_phi_clusters = tpot_container->n_phi_clusters;
    int n_phi_clusters = 0;
    
    std::vector<unsigned short> n_phi_module_clusters;
    for( int tile = 0; tile < 8; ++tile ) 
    { 
      n_phi_module_clusters.push_back( tpot_container->n_detector_clusters[tile] ); 
      if( tile > 0 ) n_phi_clusters+=tpot_container->n_detector_clusters[tile];
    }
    
    tpot_tree->GetEntry(tpot_entry+ioffset);
    // const auto n_z_clusters = tpot_container->n_z_clusters;
    int n_z_clusters = 0;

    std::vector<unsigned short> n_z_module_clusters;
    for( int tile = 0; tile < 8; ++tile ) 
    {
      n_z_module_clusters.push_back( tpot_container->n_detector_clusters[tile+8] ); 
      n_z_clusters+=tpot_container->n_detector_clusters[tile+8];
    }

    h_correlation->Fill( n_phi_clusters, n_z_clusters );

    for( int tile = 0; tile < 8; ++tile ) 
    { h_module_correlation[tile]->Fill( n_phi_module_clusters[tile], n_z_module_clusters[tile] ); }
    
  }
    
  if( true )
  {
    TCanvas* cv = new TCanvas( "cv", "cv", 1200, 900 );
    h_correlation->Draw("colz");
    // h_correlation->GetXaxis()->SetNdivisions(505);
    gPad->SetLeftMargin( 0.14);
    gPad->SetRightMargin( 0.13);
    gPad->SetLogz(true);
    
    
    auto text = new TPaveText(0.57,0.16,0.86,0.36, "NDC" );
    text->SetFillColor(0);
    text->SetFillStyle(1010);
    text->SetBorderSize(0);
    text->SetTextAlign(11);
    text->AddText("" );
    text->AddText("#it{#bf{sPHENIX}} Internal");
    text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
    text->AddText(Form("07/11/2023, Run %i", runnumber ));
    text->Draw();
    pdfDocument.Add(cv);
    cv->SaveAs(pngfilename);
  }

  if( true )
  {
    TCanvas* cv = new TCanvas( "cv1", "cv", 1200, 600 );
    cv->Divide( 4, 2 );
    for( int tile = 0; tile < 8; ++tile ) 
    { 
      cv->cd( tile+1 );
      h_module_correlation[tile]->Draw("colz");
      gPad->SetRightMargin( 0.15);
      gPad->SetLogz(true);
    }
    pdfDocument.Add(cv);
  }

}

//_______________________________________________________-
void process_all()
{
  for( const int runnumber: {20445, 20446} )
  { TPOT_Correlation( runnumber ); }
}
