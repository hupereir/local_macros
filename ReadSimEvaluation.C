
R__LOAD_LIBRARY(libg4eval_hp_lite.so)


namespace
{

  // global variables
  SimEvaluator_hp::Container* container = nullptr;
  TTree* tree = nullptr;

  // square
  double square(double x)
  { return x*x; }

  // radius
  double get_r(double x, double y)
  { return std::sqrt(square(x)+square(y)); }
}

//____________________________________________________________________________
void DrawEvent( int i )
{

  // load first event
  tree->GetEntry(i);

  // create canvas
  TCanvas* cv = new TCanvas( "cv", "", 900, 450 );
  cv->Divide( 2, 1 );
  cv->SetTopMargin(0.05);

  {
    // x,y view
    cv->cd(1);

    // create empty histogram to show axis
    TH2* h0 = new TH2I( "h0", "", 200, -100, 100, 200, -100, 100 );
    h0->GetXaxis()->SetTitle( "x (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->SetStats(0);
    h0->Draw();

    // create generic marker to show clusters
    TMarker* marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(0.6);
    marker->SetMarkerColor(2);

    // loop over container tracks
    for( const SimEvaluator_hp::G4HitStruct& g4hit:container->g4hits() )
    {
      // skip calorimeter
      if( g4hit._caloid >= 0 ) { continue; }
      marker->DrawMarker( g4hit._x, g4hit._y );
    }

  }

  {
    // r,z view
    cv->cd(2);

    // create empty histogram to show axis
    TH2* h0 = new TH2I( "h1", "", 200, -100, 100, 200, -100, 100 );
    h0->GetXaxis()->SetTitle( "z (cm)" );
    h0->GetYaxis()->SetTitle( "y (cm)" );
    h0->SetStats(0);
    h0->Draw();

    // create generic marker to show clusters
    TMarker* marker( new TMarker() );
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(0.6);
    marker->SetMarkerColor(2);

    // loop over container tracks
    for( const SimEvaluator_hp::G4HitStruct& g4hit:container->g4hits() )
    {
      // skip calorimeter
      if( g4hit._caloid >= 0 ) { continue; }

      marker->DrawMarker( g4hit._z, g4hit._y );
    }
   }

}

//____________________________________________________________________________
void ReadSimEvaluation()
{
  const TString inputFile = "DST/CONDOR_upsilon_pythia8_tau/sim_evaluation/dst_eval_sim_evaluation_upsilon_pythia8_tau_0001.root";

  auto tfile = TFile::Open( inputFile );

  // load tree and assign branch
  tree = static_cast<TTree*>( tfile->Get("T") );
  tree->SetBranchAddress( "DST#EVAL#SimEvaluator_hp::Container", &container );

  DrawEvent(0);

}
