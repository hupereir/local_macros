#include <TCanvas.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>

void DrawTPC()
{

  constexpr float rMin = 20;
  constexpr float rMax = 80;
  constexpr float length = 210;
  constexpr int nSlices = 50;

  auto manager = new TGeoManager( "tube", "poza2");
  auto material = new TGeoMaterial( "Al", 26.98,13,2.7);
  auto medium = new TGeoMedium( "MED",1,material);
  auto top = manager->MakeBox( "TOP", medium, 150, 150, 150);
  manager->SetTopVolume(top);

  auto tube = gGeoManager->MakeTube( "TUBE", medium, rMin, rMax, length/2/nSlices );
  tube->SetLineWidth(1);

  TGeoRotation rotation;
  rotation.SetAngles( 30, -140, 0 );

  for( Int_t i = 0; i < nSlices; ++i )
  {
    const float offset = -length + i*length/nSlices + length/nSlices/2;
    TGeoTranslation translation( 0, 0, offset );
    top->AddNode( tube, i+1, new TGeoHMatrix( rotation*translation ) );
  }

  manager->CloseGeometry();
  manager->SetNsegments(72);

  auto cv = new TCanvas("c", "c",0,0,600,600);
  top->Draw();

}
