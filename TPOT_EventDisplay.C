#include <TVector3.h>

#include <optional>

#include <nlohmann/json.hpp>

// TPOT cluster
class Cluster
{
  public:

  unsigned short layer = 0;
  
  unsigned short tile = 0;

  unsigned short size = 0;
  
  double charge = 0;

  int strip = 0;

  int region = 0;
  
  TVector3 begin;
  TVector3 end;  
  TVector3 center;  

  bool is_signal = true;
  
  using List = std::vector<Cluster>;
};

// TPOT event
class FullEventContainer: public TObject
{
  
  public:
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
  
  // clusters
  Cluster::List clusters;
  
  // number of clusters per detector
  std::vector<unsigned short> n_detector_clusters;
 
  // number of clusters per region
  std::vector<unsigned short> n_region_clusters;
 
  // minimum cluster charge per detector
  std::vector<unsigned short> min_cluster_size;
 
  // minimum cluster charge per detector
  std::vector<double> min_cluster_charge;
  
  // number of clusters
  unsigned short n_phi_clusters = 0;
  
  // number of clusters
  unsigned short n_z_clusters = 0;
    
  ClassDef(FullEventContainer,1)
  
};

namespace
{
  // small class to easily save multple pages pdf
  class PdfDocument
  {
    public:
    
    //* constructor
    PdfDocument( TString filename ): m_filename( filename ) {}    
    
    //* destructor
    ~PdfDocument()
    {
      if( m_first || m_filename.IsNull() ) return;
      // we save a new empty canvas at the end of the pdf file so that it is properly closed
      TCanvas().SaveAs( Form( "%s)", m_filename.Data() ) );
    }
    
    //* add pad
    void Add( TVirtualPad* pad )
    {
      if( m_filename.IsNull() ) return;
      if( m_first )
      {        
        pad->SaveAs( Form( "%s(", m_filename.Data() ) );
        m_first = false;
      } else {
        pad->SaveAs( Form( "%s", m_filename.Data() ) );
      }
    }
    
    private:
    
    //* filename
    TString m_filename;
    
    //* true if first page
    Bool_t m_first = true;
  };
  
}

namespace 
{
  using line_t = std::array<double, 4 >;
  
  //____________________________________________________________________________
  void draw_tpot_xy()
  {
    
    std::vector<line_t> tpot_xy_lines = 
    { 
      { -14.4999, -83.3787, 10.9999, -83.4836}, 
      { -14.4849, -83.2659, 11.0151, -83.3088}, 
      { -14.4153, -83.0462, 11.0828, -83.354}, 
      { -14.4052, -83.3156, 11.0945, -83.4362}, 
      { -54.6824, -65.0237, -32.4596, -77.5297}, 
      { -54.6245, -64.9096, -32.425, -77.4569}, 
      { 29.2854, -79.0461, 51.2348, -66.0665}, 
      { 29.4241, -79.0645, 51.3265, -66.0058}
    }; 

    auto line = new TLine();
    line->SetLineWidth(2);
    line->SetLineColor(1);
    
    for( const auto& [x1, y1, x2, y2]:tpot_xy_lines )
    { line->DrawLine( x1, y1, x2, y2 ); }
  }

  //____________________________________________________________________________
  void draw_tpot_rz()
  {
        
    std::vector<line_t> tpot_rz_lines = { 
      { -110.041, -84.0957, -59.0417, -83.9141}, 
      { -53.5958, -84.0698, -2.59725, -83.6961}, 
      { 2.7824, -83.7376, 53.7822, -83.8303}, 
      { 59.1906, -84.0373, 110.191, -84.0589}, 
      { -62.2889, -84.1996, -11.2893, -83.9958}, 
      { 11.8928, -83.9804, 62.8928, -83.9946}, 
      { -62.8673, -83.7926, -11.8677, -83.6068}, 
      { 11.3459, -83.6263, 62.3457, -83.6154}
    }; 
    
    auto line = new TLine();
    line->SetLineWidth(2);
    line->SetLineColor(1);
    
    for( const auto& [z1, r1, z2, r2]:tpot_rz_lines )
    { line->DrawLine( z1, r1, z2, r2 ); }
  }

}

namespace
{
  
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
    
  template<class T>
    inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  // global run number
  int m_runnumber = 0;
  
  FullEventContainer* tpot_event = nullptr;
    
  //_______________________________________________
  void setup_tpot_tree( TTree* tree )
  { tree->SetBranchAddress( "Event", &tpot_event ); }  
  
  //____________________________________________________________________________
  void write_to_json( const TString& filename )
  {

    nlohmann::json runstats;
    runstats.push_back( "sPHENIX Experiment at RHIC" );
    runstats.push_back( "Data recorded: 2023-08-11" );
    runstats.push_back( Form( "Run / Event: %i/%i", m_runnumber,  tpot_event->lvl1_counter ) );
    runstats.push_back( "Cosmics" );
    nlohmann::json event;
    event["runstats"] = runstats;
    event["B"] = 1.4;
    
    // meta data
    nlohmann::json meta = 
    { { "HITS",{
        { "INNERTRACKER", 
          {
            { "type", "3D" },
            { "options", { {"size", 10}, { "color", 16777215 } } }
          }
        },
        { "TRACKHITS",
          {
            { "type", "3D" },
            { "options", { {"size", 2 }, { "transparent", 0.8 }, { "color", 16777215 } } }
          }
        }
      } }
    };
        
    nlohmann::json hits;    
    
    // add TPOT clusters
    static constexpr int npoints = 100;
    for( const auto& cluster:tpot_event->clusters )
    {
      for( int ipoint = 0; ipoint <= npoints; ++ipoint )
      {
        double fraction = double(ipoint)/npoints;
        const auto strip_center = fraction*cluster.begin + (1.0-fraction)*cluster.end;
     
        nlohmann::json strip_info = {
          {"x", strip_center.x()},
          {"y", strip_center.y()},
          {"z", strip_center.z()},
          {"e", 1 }
        };
        hits["TRACKHITS"].emplace_back( strip_info );
      }
      
    }
   
    nlohmann::json main;
    main["EVENT"]=event;
    main["META"]=meta;
    main["HITS"]=hits;            
  
    // write event to file
    std::ofstream out( filename.Data() );
    out << std::setw(2) << main << std::endl;
    out.close();  
        
  }
}

//____________________________________________________________________________
TString TPOT_EventDisplay( int runnumber = 25926 )
{
  // store run number globally. Needed for JSon
  m_runnumber = runnumber;
 

  // pdf output
  const TString pdffilename = Form( "Figures/TPOT_event_display-%08i-0000.pdf", runnumber );
  PdfDocument pdfDocument( pdffilename );
  
  // tpot input file
  const TString tpot_inputfilename = Form( "Rootfiles/RawDataClusterTree-%08i-0000-test_filtered.root", runnumber );
  auto tpot_tfile = TFile::Open( tpot_inputfilename, "READ" );
  auto tpot_tree = static_cast<TTree*>( tpot_tfile->Get("T") );
  setup_tpot_tree( tpot_tree );
  std::cout << "TPOT_EventDisplay - tpot_inputfilename: " << tpot_inputfilename << std::endl;
  std::cout << "TPOT_EventDisplay - pdffilename: " << pdffilename << std::endl;

  const auto tpot_entries = tpot_tree->GetEntries();
  int tpot_entry = 0;
  for(; tpot_entry < tpot_entries; ++tpot_entry )
  {

    tpot_tree->GetEntry( tpot_entry );
    auto tpot_bco = tpot_event->lvl1_bco;
    
    // create canvas
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    cv->cd(1);

    // std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -15, 15, 100, -15, 15 ) );
    std::unique_ptr<TH2> h0( new TH2I( "h0", "", 100, -60, 60, 100, -95, 15 ) );
    h0->SetTitle( Form( "runnumber: %i, event: %i, BCO: %lu", runnumber, tpot_event->lvl1_counter, tpot_bco ) );
    h0->GetXaxis()->SetTitle( "x (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->Draw();
    
    std::unique_ptr<TMarker> marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerColor(2);
        
    // draw TPOT clusters
    draw_tpot_xy();
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 55 )
      { marker->DrawMarker( cluster.center.x(), cluster.center.y() ); }
    }
        
    cv->cd(2);
    // std::unique_ptr<TH2> h1( new TH2I( "h1", "", 100, -50, 50, 100, -15, 15 ) );
    std::unique_ptr<TH2> h1( new TH2I( "h1", "", 100, -110, 110, 100, -95, 15 ) );
    h1->GetXaxis()->SetTitle( "z (cm)" );
    h1->GetYaxis()->SetTitle( "r (cm)" );
    h1->Draw();
              
    // draw TPOT clusters  
    draw_tpot_rz();
    for( const auto& cluster:tpot_event->clusters )
    {
      if( cluster.layer == 56 )
      { marker->DrawMarker( cluster.center.z(), -get_r(cluster.center.x(), cluster.center.y()) ); }
    }

    pdfDocument.Add( cv.get() );
    
    // generate json file
    const auto json_filename = Form( "event_display/tpot_event_display_%08i_%08i_%lu.json", runnumber, tpot_event->lvl1_counter, tpot_bco );
    write_to_json( json_filename);
  }
  
  return pdffilename;
  
}
