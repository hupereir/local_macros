
#include <TCanvas.h>

void DrawEfficiency_Saclay()
{

    const int n_points_phi  = 5;
    double resist_phi[5] = { 420, 410, 400, 380, 350 };
    double eff_phi[5] = { 0.99516, 0.99403, 0.99014, 0.9451, 0.62991 };
    auto tg_phi = new TGraph( n_points_phi, resist_phi, eff_phi );
    tg_phi->SetMarkerStyle(20);
    tg_phi->SetMarkerColor(kBlue);
    tg_phi->SetLineColor(kBlue);
    tg_phi->SetMarkerSize(2);

    const int n_points_z  = 5;
    const double resist_z[5] = { 350, 380, 400, 410, 420 };
    const double eff_z[5] = { 0.51869, 0.89326, 0.96707, 0.98046, 0.97807 };
    auto tg_z = new TGraph( n_points_z, resist_z, eff_z );
    tg_z->SetMarkerStyle(20);
    tg_z->SetMarkerColor(kRed);
    tg_z->SetLineColor(kRed);
    tg_z->SetMarkerSize(2);


    // create plot
    TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );

    auto h = new TH1I( "h", "", 100, 310, 440 );
    h->SetMinimum(0);
    h->SetMaximum(1);
    h->GetXaxis()->SetTitle( "resist HV (V)" );
    h->GetYaxis()->SetTitle( "efficiency" );
    h->Draw();
    gPad->SetTopMargin( 0.05 );
    gPad->SetLeftMargin( 0.14);

    gPad->Update();
    tg_phi->Draw("P");
    tg_z->Draw("P");

}
