//___________________________________________________________
Double_t Function( Double_t* x, Double_t* par)
{
  Double_t alpha = (x[0]-par[1])/par[2];
  Double_t alpha2 = -(x[0]-par[6])/par[7];
  return
    par[0]*(
    par[3]*(TMath::Erf( alpha/TMath::Sqrt( 2 ) )+1)/2 +
    (1.0-par[3])*(TMath::Erf( alpha/par[4]/TMath::Sqrt( 2 ) )+1)/2)
    +par[5]*(TMath::Erf( alpha2/TMath::Sqrt( 2 ) )+1)/2
    ;
}

//___________________________________________________________
Double_t Function( Double_t* x, Double_t* par )
{

  Double_t alpha1 = (x[0]-par[1])/par[2];
  Double_t alpha2 = -(x[0]-par[4])/par[5];
  return
    par[0]*(TMath::Erf( alpha1/TMath::Sqrt( 2 ) )+1)/2 +
    par[3]*(TMath::Erf( alpha2/TMath::Sqrt( 2 ) )+1)/2;

}

void DrawErf( void )
{
  // TF1* f = new TF1( "Function", Function, -5, 5, 5 );
  // f->SetParameters( 1, 0.5, 0.3, 1.0, 2 );

  // TF1* f = new TF1( "Function", Function, -5, 5, 8 );
  // f->SetParameters( 1, 0.5, 0.3, 1.0, 2, 1, 0, 0.5 );

  TF1* f = new TF1( "Function", Function, -5, 5, 6 );
  f->SetParameters( 1.0, 1., 0.3, 1., 0.0, 0.1);

  f->SetLineColor(2);
  f->Draw();
}
