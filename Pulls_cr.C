
#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>

//____________________________________________________________________________
Float_t delta_phi( Float_t phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI; 
  else return phi;
}

//____________________________________________________________________________
void Pulls_cr()
{
  
  const TString inputFile = "dst_eval_cr.root";
  const TString pdfFile = "Pulls_cr.pdf";
  std::cout << "Pulls_cr - inputFile: " << inputFile << std::endl;
  std::cout << "Pulls_cr - pdfFile: " << pdfFile << std::endl;
  
  auto tfile = TFile::Open( inputFile );
  if( !tfile ) {
    std::cout << "Pulls_cr - invalid input file: " << inputFile << std::endl;
    return;
  }
  
  auto tree = dynamic_cast<TTree*>( tfile->Get("T") );
  if( !tree )
  {
    std::cout << "Pulls_cr - invalid tre \"T\"" << std::endl;
    return;
  }
  
  // define variables and cuts
  /* this is cluster - truth pulls, for clusters on tracks */
  const TString var( "delta_phi(m_tracks.clusters.phi - m_tracks.clusters.truth_phi)/m_tracks.clusters.phi_error" );

  /* same thing with z */
  // const TString var( "(m_tracks.clusters.z - clusters.truth_z)/m_tracks.clusters.z_error" );
  
  /* 
   * this is cluster - track pulls, for clusters on tracks  
   * in principle one should also add the track error quadratically to the denominator.
   * It is expected to be small and its validity should be checked anyway
   */
  // const TString var( "delta_phi(m_tracks.clusters.phi - m_tracks.clusters.trk_phi)/m_tracks.clusters.phi_error" );

  /* same thing with z */  
  // const TString var( "m_tracks.clusters.z - m_tracks.clusters.trk_z)/m_tracks.clusters.z_error" );

  // layer selection
  /* 23 is the first layer of the R2 GEMS in TPC */
  const int ilayer = 23;
  const TCut layer_cut = Form( "m_tracks.clusters.layer == %i", ilayer );
  
  // cluster size cut
  /* no cut */
  const TCut csize_cut = ("m_tracks.clusters.phi_size > 0 && m_tracks.clusters.z_size > 0" );
  
  // create histogram
  const TString hname = "h_pullsrph";
  auto h = new TH1F( hname, "", 100, -5, 5 );
  h->GetXaxis()->SetTitle( "#Delta#phi_{clus-truth}/#sigma#phi" );
  // h->GetXaxis()->SetTitle( "#Deltaz_{clus-truth}/#sigmaz" );
  
  // project tree
  tree->Project( hname, var, layer_cut && csize_cut );

  // fit
  h->Fit( "gaus", "0" );
  
  // draw
  auto cv = new TCanvas( "cv", "cv", 900, 900 );
  h->Draw();
  
  auto f = h->GetFunction("gaus");
  f->SetLineColor(2);
  f->Draw("same");
  
  cv->SaveAs( pdfFile );
} 
