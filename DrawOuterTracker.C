#include <TCanvas.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>

namespace 
{
  constexpr float rMin = 20;
  constexpr float rMax = 80;
  constexpr float length = 210;
  
  constexpr float mm_radius = 82;
  constexpr float mm_thickness = 2;
  constexpr float mm_length = 50;
  constexpr float mm_width = 25;
  
  TGeoManager* manager = nullptr;
  TGeoMedium* medium = nullptr;
  TGeoVolume* rotated = nullptr;
}


// TPC volume
void DrawTPC()
{
  auto tube = manager->MakeTube( "TUBE", medium, rMin, rMax, length/2 );
  tube->SetLineWidth(1);
  tube->SetLineColor( 18 );
  rotated->AddNode( tube, 0 );
}

// TPC readout
void DrawTPCReadout()
{
  auto gem = manager->MakeTubs( "GEM", medium, rMin, rMax-1, 0.5, -15, 15 );
  gem->SetLineColor( kCyan-3);
    
  for( int i = 0; i < 12; ++i )
  {
    TGeoTranslation translation( 0, 0, length/2 );
    TGeoRotation rotation;
    rotation.RotateZ(15+30*i);
    rotated->AddNode( gem, i+1, new TGeoHMatrix(rotation*translation) );
  }
  
  
  for( int i = 0; i < 12; ++i )
  {
    TGeoTranslation translation( 0, 0, -length/2 );
    TGeoRotation rotation;
    rotation.RotateZ(15+30*i);
    rotated->AddNode( gem, i+1, new TGeoHMatrix(rotation*translation) );
  }
  
}

// TPC strongbackreadout
void DrawStrongBack()
{
  constexpr double half_thickness = 3;
  constexpr float rGemMax = 79.5;
  auto separator = manager->MakeTrd1( "SEP", medium, 1, 1, half_thickness, (rGemMax-rMin)/2 );
  separator->SetLineColor(18);
  
  TGeoRotation rotation1;
  rotation1.RotateX(-90);
  
  for( int i = 0; i < 12; ++i )
  {
    TGeoTranslation translation( 0, (rGemMax+rMin)/2, length/2+half_thickness );
    TGeoRotation rotation;
    rotation.RotateZ(30*i);
    rotated->AddNode( separator, i+1, new TGeoHMatrix(rotation*translation*rotation1) );
  }
  
  
  for( int i = 0; i < 12; ++i )
  {
    TGeoTranslation translation( 0, (rGemMax+rMin)/2, -length/2-half_thickness );
    TGeoRotation rotation;
    rotation.RotateZ(30*i);
    rotated->AddNode( separator, i+1, new TGeoHMatrix(rotation*translation*rotation1) );
  }
  
  static constexpr std::array<double,4> radius = {{20, 40, 60, 80}};
  for( int i = 0; i < 4; ++i )
  {
    auto tube = manager->MakeTube( Form("TUBE%i", i ), medium, radius[i]-0.5, radius[i]+0.5, half_thickness );
    tube->SetLineColor(18);

    {
      TGeoTranslation translation( 0, 0, +length/2+half_thickness );
      rotated->AddNode( tube, 1, new TGeoHMatrix(translation) );
    }
    
    {
      TGeoTranslation translation( 0, 0, -length/2-half_thickness );
      rotated->AddNode( tube, 2, new TGeoHMatrix(translation) );
    }
    
  }
  
}

TString DrawMicromegas0()
{
  const auto angle_offset = 15;
  const auto angle = (180./M_PI)*(mm_width/mm_radius);
  auto box = manager->MakeTubs( "BOX", medium, mm_radius, mm_radius+mm_thickness, mm_length/2, -angle/2 + angle_offset, angle/2 + angle_offset );
  box->SetLineColor( kCyan-3);

  {
    TGeoTranslation translation( 0, 0, mm_length/2 + 10 );
    rotated->AddNode( box, 1, new TGeoHMatrix( translation ) ); 
  }
  
  {
    TGeoTranslation translation( 0, 0, -mm_length/2 - 10 );
    rotated->AddNode( box, 2, new TGeoHMatrix( translation ) ); 
  }
  return "outertracker-0.png";
}

// micromegas four thin sector (curved)
TString DrawMicromegas1()
{
  const auto angle_offset = 15;
  const auto angle = (180./M_PI)*(mm_width/mm_radius);
  auto box = manager->MakeTubs( "BOX", medium, mm_radius, mm_radius+mm_thickness, mm_length/2, -angle/2 + angle_offset, angle/2 + angle_offset );
  box->SetLineColor( kCyan-3);
  for( int i = 0; i < 4; ++i )
  {
    TGeoTranslation translation( 0, 0, -length/2+(i*length/4) + mm_length/2 );
    rotated->AddNode( box, i+1, new TGeoHMatrix( translation ) );
  }
  
  return "outertracker-1.png";

}   

// micromegas 12 thin sectors (curved)
TString DrawMicromegas2()
{
  const auto angle_offset = 15;
  const auto angle = (180./M_PI)*(mm_width/mm_radius);
  auto box = manager->MakeTubs( "BOX", medium, mm_radius, mm_radius+mm_thickness, mm_length/2, -angle/2 + angle_offset, angle/2 + angle_offset );

  box->SetLineColor( kCyan-3);
  for( int i = 0; i < 12; ++i )
  {
    TGeoRotation rotation;
    rotation.RotateZ(30*i);
    
    rotated->AddNode( box, i+1, new TGeoHMatrix( rotation ) );
  }
  return "outertracker-2.png";
}

// micromegas 24 thin sectors (curved)
TString DrawMicromegas3()
{
  const auto angle_offset = 15;
  const auto angle = (180./M_PI)*(mm_width/mm_radius);
  auto box = manager->MakeTubs( "BOX", medium, mm_radius, mm_radius+mm_thickness, mm_length/2, -angle/2 + angle_offset, angle/2 + angle_offset );

  box->SetLineColor( kCyan-3);
  for( int i = 0; i < 12; ++i )
  {
    TGeoRotation rotation;
    rotation.RotateZ(30*i);
    
    {
      TGeoTranslation translation( 0, 0, mm_length/2 + 10 );
      rotated->AddNode( box, 2*i+1, new TGeoHMatrix( translation*rotation ) ); 
    }

    {
      TGeoTranslation translation( 0, 0, -mm_length/2 - 10 );
      rotated->AddNode( box, 2*i+2, new TGeoHMatrix( translation*rotation ) ); 
    }
  }
  return "outertracker-3.png";
}

void DrawOuterTracker( int selection )
{
  
  TString picturename = "outertracker.png";
  manager = new TGeoManager( "TOC", "TPC");
  
  // main volume
  auto matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  auto vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  auto top = manager->MakeBox( "TOP", vacuum, 150, 150, 150);
  manager->SetTopVolume(top);

  rotated = manager->MakeBox( "Rotated", vacuum, 150, 150, 150);
  rotated->SetVisibility( false );
  {
    TGeoRotation rotation;
    rotation.RotateZ( 120 );
    rotation.RotateX( -10 );
    rotation.RotateY( 40 );
  
    top->AddNode( rotated, 0, new TGeoHMatrix( rotation ) );
  }
  
  auto material = new TGeoMaterial( "Al", 26.98,13,2.7);
  medium = new TGeoMedium( "MED",1,material);
  
  // TPC volume
  DrawTPC();
  DrawTPCReadout();
  DrawStrongBack();
  
  // Outer Tracker
  switch( selection )
  {
    case 0: picturename = DrawMicromegas0(); break;
    case 1: picturename = DrawMicromegas1(); break;
    case 2: picturename = DrawMicromegas2(); break;
    case 3: picturename = DrawMicromegas3(); break;
    default: break;
  }
  
  manager->CloseGeometry();

  top->Draw( "ogl" );
  
  auto viewer = static_cast<TGLViewer*>(gPad->GetViewer3D());
  viewer->SavePicture( picturename );
}

void DrawOuterTracker()
{ 
  DrawOuterTracker(0); 
  DrawOuterTracker(1); 
  DrawOuterTracker(2); 
   DrawOuterTracker(3); 
}
