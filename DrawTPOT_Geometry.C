
#include <TCanvas.h>
#include <TH1.h>
#include <TLine.h>

#include <array>
#include <vector>


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


//____________________________________________________________________________
void DrawTPOT_Geometry()
{
  
  
  // create canvas
  auto cv( new TCanvas( "cv", "cv", 1200, 600 ) );
  cv->Divide( 2, 1 );
  
  
  cv->cd(1);
  auto h0( new TH2I( "h0", "", 100, -60, 60, 100, -95, 15 ) );
  h0->GetXaxis()->SetTitle( "x (cm)" );
  h0->GetYaxis()->SetTitle( "y (cm)" );
  h0->Draw();
  draw_tpot_xy();

  cv->cd(2);
  auto h1( new TH2I( "h1", "", 100, -110, 110, 100, -95, 15 ) );
  h1->GetXaxis()->SetTitle( "z (cm)" );
  h1->GetYaxis()->SetTitle( "r (cm)" );
  h1->Draw();
  draw_tpot_rz();
}
 
