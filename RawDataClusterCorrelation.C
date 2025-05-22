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
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name ) );
    }
  }

  int get_color( int layer, int tile, int region )
  {
    std::array<int, 8> color_base = {kOrange, kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan };
    if( layer == 0 ) return color_base[tile]+region;
    else return color_base[tile]-1-region;
  }


  //_______________________________________________________________
  TLine* draw_diagonal( double min, double max )
  {
    auto line = new TLine( min, min, max, max);
    line->SetLineStyle( 2 );
    line->SetLineWidth( 1 );
    line->SetLineColor( 1 );
    return line;
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
}


//__________________________________________________________________________________
void RawDataClusterCorrelation(int runnumber = 54535)
{

  const auto inputfile = Form( "DST/dst_eval-%08i-0000.root", runnumber );
  const auto pdffile = Form( "Figures/RawDataClusterCorrelation-%08i-0000.pdf", runnumber );

  std::cout << "RawDataClusterCorrelation - inputfile: " << inputfile << std::endl;
  std::cout << "RawDataClusterCorrelation - pdffile: " << pdffile << std::endl;

  const auto run_label = make_run_label({runnumber});
  const auto postfix = make_run_postfix({runnumber});

  PdfDocument pdfdocument(pdffile);

  save_detector_names();

  FileManager fileManager( inputfile );
  auto tree = fileManager.GetChain( "T" );

//   Utils::max_entries = 10000;
  if( true )
  {
    // per tile correlation between phi vs z
    static const int max_clusters = 30;

    auto cv = new TCanvas( "cv0", "cv0", 1200, 620 );
    cv->Divide( 4, 2, 0.002, 0.002 );

    // one canvas per tile
    for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile;
      const int detid_other = itile+8;

      const auto& detname = detector_names[detid];

      const auto hname = Form( "h_tile_%i", detid );
      auto h = new TH2I( hname, "", max_clusters-1, 1, max_clusters, max_clusters-1, 1, max_clusters );
      auto var = Form( "n_detector_clusters[%i]:n_detector_clusters[%i]", detid_other, detid );
      Utils::TreeToHisto( tree, hname, var, TCut(), false );

      cv->cd( detid+1 );
      h->SetStats(0);
      h->SetTitle("");
      h->GetXaxis()->SetTitle(Form("N_{clus} %s",  detector_names[detid].c_str()));
      h->GetYaxis()->SetTitle(Form("N_{clus} %s",  detector_names[detid_other].c_str()));
      h->Draw("colz");
      gPad->SetLogz(true);

      gPad->SetTopMargin(0.1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.13);
      Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.substr(0,3).c_str(), run_label.Data(), h->GetEntries() ));



      draw_diagonal(1,max_clusters)->Draw();
    }

    pdfdocument.Add( cv );
  }

  if( true )
  {

    // correlation between all phi vs all z
    auto cv = new TCanvas( "cv1", "cv1", 800, 800 );

    static const int max_clusters = 180;
    const TString hname( "h_all" );
    auto h = new TH2I( hname, "", max_clusters-1, 1, max_clusters, max_clusters-1, 1, max_clusters );

    const TString var = "(n_detector_clusters[9]+n_detector_clusters[10]+n_detector_clusters[11]+n_detector_clusters[12]+n_detector_clusters[13]+n_detector_clusters[14]+n_detector_clusters[15]):"
      "(n_detector_clusters[1]+n_detector_clusters[2]+n_detector_clusters[3]+n_detector_clusters[4]+n_detector_clusters[5]+n_detector_clusters[6]+n_detector_clusters[7])";

    Utils::TreeToHisto( tree, hname, var, TCut(), false );
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("N_{clus} (#phi views)" );
    h->GetYaxis()->SetTitle("N_{clus} (z views)" );
    h->GetYaxis()->SetTitleOffset(1.4);
    h->Draw("colz");
    gPad->SetLogz(true);

    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.13);
    Draw::PutText( 0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h->GetEntries() ));

    draw_diagonal(1,max_clusters)->Draw();

    pdfdocument.Add( cv );
  }

  if( true )
  {

    // correlation between all phi vs all z
    auto cv = new TCanvas( "cv2", "cv2", 1000, 800 );

    static const int max_charge = 2000;
    static const int max_clusters = 360;
    const TString hname( "h_all" );
    auto h = new TH2I( hname, "", 100, 0, max_charge, max_clusters, 0, max_clusters );

    const TString var = "n_clusters:mbd_info.charge_south+mbd_info.charge_north";

    Utils::TreeToHisto( tree, hname, var, TCut(), false );
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("MBD Q_{tot}");
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle("TPOT N_{clusters}" );
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw("colz");
    gPad->SetLogz(true);

    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.125);
    gPad->SetRightMargin(0.125);
    Draw::PutText( 0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h->GetEntries() ));

    pdfdocument.Add( cv );
  }
}
