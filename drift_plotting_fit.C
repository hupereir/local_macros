#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <TChain.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TFile.h>

void append_to_csv(std::string filename, std::string data)
{
  std::ifstream infile(filename);
  bool fileExists = infile.good();
  infile.close();
  std::ofstream file(filename, std::ios::app);
  if (file.is_open())
    {
      if(!fileExists)
	{
	  file << "Run Number,Number of Entries Processed,Processed Drift Velocity (cm/ns),New Drift Velocity (cm/ns),Drift Velocity Error (cm/ns),New T0 Correction(bins),T0 Correction Error(bins),z Misalignment (cm),z Misalignment Error (cm),Total Jobs Processed (10 Jobs Per Segment),Zombie Jobs"<<"\n";
	}
      file << data;
      file << "\n";
      file.close();
    }
  else
    {
      std::cerr << "Can't open the file " << filename << std::endl;
    }
}

std::map<int, float> csv_to_map(std::string filename, int column)
{
  std::ifstream file(filename);
  std::map<int, float> data;
  std::string line;

  //skip the first line
  std::getline(file, line);
  while (std::getline(file, line))
    {
      std::stringstream ss(line);
      std::string token;
      int runNumber;
      float value;
      int icol = 0;

      while (std::getline(ss, token, ','))
	{
	  if (icol == 0) runNumber = std::stoi(token);
	  if (icol == column) value = std::stof(token);
	  icol++;
	}

      data[runNumber] = value;
    }

  file.close();
  return data;

}

std::map<int, std::string>get_module_names()
{
  std::map<int, std::string>modules={{0,"SCOZ"},{1,"SCIZ"},{2,"NCIZ"},{3,"NCOZ"},{4,"SEZ"}, {5,"NEZ"}, {6,"SWZ"}, {7,"NWZ"}};
  return modules;
}

std::vector<std::string> read_filenames(const char* filename)
{
  std::ifstream file(filename);
  std::vector<std::string> filenames;

  std::string line;
  while (std::getline(file, line))
    {
      filenames.push_back(line);
    }
  file.close();
  return filenames;
}
double fit_function(double* x, double* par)
{
  const int itile = std::floor(x[0]); // Not entirely sure why we have to floor this
  const double z = x[1];

  if( itile >= 8 ) { TF1::RejectPoint(); return 0.0; }

  const double slope = par[0];
  const double offset = par[itile+1];
  return offset+slope*z;
}

double linear_function( double* x, double* par )
{
  return par[0]*x[0]+par[1];
}

void drift_plotting_step(int runnumber, float v_drift, float t0)
{

  std::vector<std::string> filePaths = read_filenames(Form("/gpfs/mnt/gpfs02/sphenix/user/bade/sphenix/work/g4eval/file_lists/file_list_%i_3rd_order_24022025.txt",runnumber));
  std::map<int, std::string>module_names = get_module_names();
  TChain *tree = new TChain("T");
  int nZombies = 0;
  int nTotal = 0;
  for (const auto& file : filePaths)
    {

      //The Zombie counter. Although it's not too necessary for this macro
      TFile* f = TFile::Open(file.c_str());
      nTotal++;
      if (!f || f->IsZombie())
	{
	  nZombies++;
	  if (f) f->Close();
	  continue;
	}

      tree->Add(file.c_str());
      f->Close();
    }
  int nEntries = tree->GetEntries();
  std::cout<<"For Run number "<<runnumber<<" nEntries = "<<nEntries<<" and number of zombies = "<<nZombies<<std::endl;

  const TCut track_cut(
		       " 1 "
		       //    "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._pt>0.2"
		       // "&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing==1 || DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing==0)"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing<1000"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._ndf > 0"
		       "&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._chisquare/DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._ndf)<100"
		       //"&& (DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._chisquare/DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._ndf)<1000"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_tpc>20"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_mvtx>2"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._nclusters_intt>1"
		       "&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._layer==56 && DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._found_cluster_z._layer==56"
		       //"&& DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._tile<4"
		       );

  std::cout<<"MM0"<<std::endl;

  //Keep the variable order consistent
  const TString var =  "DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._z-DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._found_cluster_z._z";
  TString var2d = Form("%s:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._z",var.Data());
  TString var3d = Form("%s:DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._trk_state_z._tile",var2d.Data());

  std::cout<<"Variables set"<<std::endl;

  static constexpr int nzbins = 50;
  TH3F *hist3D = new TH3F("hist3D", "; z (track) (cm); #Deltaz (track-cluster) (cm); tile ", 8, 0, 8, nzbins, -110, 110, 100, -10, 10);
  std::cout<<hist3D->GetEntries()<<std::endl; //This returns 0
  //tree->Scan(var3d.Data(), track_cut); //Looks fine.

  // tree->Draw(Form("%s>>hist3D", var3d.Data()), track_cut, "goff"); //This is the part that doesnt work. ROOT is not filling the existing histogram i have for some reason...? I know I'm not overfilling the histogram since this same thing worked before in 2D
  tree->Draw(Form("%s>>hist3D", var3d.Data()), track_cut, "goff", 500000); //This is the part that doesnt work. ROOT is not filling the existing histogram i have for some reason...? I know I'm not overfilling the histogram since this same thing worked before in 2D
  std::cout<<hist3D->GetEntries()<<std::endl; //This returns garbage.
  std::cout<<"MM1"<<std::endl;
  auto h_fit = new TH2F( "h_fit", "", 8, 0, 8, nzbins, -110, 110);

  for(int j=0;j<8;j++)
    {
      hist3D->GetXaxis()->SetRange(j+1, j+1);
      TH2* h2d = (TH2F*) hist3D->Project3D("zy");
      h2d->SetName(Form("h_%s", module_names[j].c_str()));
      h2d->SetTitle(Form("%s;z (track) (cm); #Deltaz (track-cluster) (cm)",module_names[j].c_str()));
      h2d->FitSlicesY(nullptr, 0, -1, 20);
      TH1F* h_mean = (TH1F*)gDirectory->Get(Form( "%s_1", h2d->GetName()));
      TH1F* h_rms = (TH1F*) gDirectory->Get(Form( "%s_2", h2d ->GetName()));

      for( int i = 0; i<h_mean->GetNbinsX(); ++i )
	{
	  const auto entries = h2d->Integral(i+1,i+1,1,hist3D->GetNbinsY());
	  if( entries > 0 )
	    {
	      h_fit->SetBinContent(j+1, i+1, h_mean->GetBinContent(i+1) );
	      h_fit->SetBinError(j+1, i+1, h_rms->GetBinContent(i+1)/std::sqrt(entries));
	    }
	}
    }
  TF2* fit = new TF2( Form( "fit_function_%i", runnumber), fit_function, 0, 8, -110, 110, 9);

  for(int i=0; i<9; i++)
    {
      fit->SetParameters(i, 0.0);
    }

  std::cout<<"MM2"<<std::endl;

  // Fit the function to the 3d histogram
  h_fit->Fit(fit, "0R"); //Okay so this is not the main problem

  std::cout<<"MM3"<<std::endl;

  int nFit = hist3D->GetEntries();

  double slope = fit->GetParameter(0);
  double new_drift = v_drift/(1+slope);
  double slope_err = fit->GetParError(0);
  double drift_err = v_drift / pow(1 + slope, 2) * slope_err;
  TCanvas *c = new TCanvas("c", "3dplot", 2000, 1000);
  c->Divide(4,2);

  std::cout<<"MM4"<<std::endl;
  for(int j=0;j<8;j++)
    {
      c->cd(j+1);

      hist3D->GetXaxis()->SetRange(j+1, j+1);
      TH2* h2d = (TH2F*) hist3D->Project3D("zy");
      h2d->SetName(Form("h_%s", module_names[j].c_str()));
      h2d->SetTitle(Form("%s;z (track) (cm); #Deltaz (track-cluster) (cm)",module_names[j].c_str()));
      h2d->Draw("COLZ" );
      h2d->SetStats(0);

      //Project the gaussian fit plot to be 1D
      auto h_fit_proj = h_fit->ProjectionY( Form( "%s%i", h_fit->GetName(), j ), j+1, j+1);
      h_fit_proj->SetMarkerStyle(20);
      h_fit_proj->SetMarkerColor(kRed);
      h_fit_proj->SetLineColor(1);
      h_fit_proj->Draw("same");

      auto f1d = new TF1( Form( "fit_function_%i", j), linear_function, -110, 110, 2 );
      f1d->SetParameter(0, slope);
      f1d->SetParameter(1, fit->GetParameter(j+1));
      f1d->SetLineColor(kGreen);
      f1d->SetLineWidth(2);
      f1d->Draw("same");

      // Create legend
      TLegend *leg = new TLegend(0.4, 0.75, 0.9, 0.9);
      leg->SetHeader( Form("%i Entries, Drift Velocity = %.2f m/ms",nFit , v_drift*100000), "C");
      leg->AddEntry(h_fit_proj,"Gaussian Fit", "p");
      leg->AddEntry(f1d, Form("Slope = %.4f", slope), "l");
      leg->AddEntry(f1d, Form("New drift velocity = %.2f#pm%.2f", new_drift*100000, drift_err*100000));
      leg->Draw();
    }

    // c->SaveAs(Form("/gpfs/mnt/gpfs02/sphenix/user/bade/sphenix/work/g4eval/plots/slicefit_pertile_t0%+.2f_%i_02032025.png", t0, runnumber));
  c->SaveAs(Form("Figures/slicefit_pertile_t0%+.2f_%i_02032025.png", t0, runnumber));

  if(false)
    {
      TCanvas *x = new TCanvas("x","crossing", 1000, 1000);
      TH1F *crossing = new TH1F("crossinghist", Form("Run %i with Silicon Match, T0%+.2f, Corrected Drift, SAMPA bias = 0ns;Crossing; Count", runnumber, t0 ), 500, -100, 400);
      tree->Draw("DST#EVAL#MicromegasTrackEvaluator_hp::Container._tracks._crossing>>crossinghist", track_cut);
      gPad->SetLogy();
      x->SaveAs(Form("/gpfs/mnt/gpfs02/sphenix/user/bade/sphenix/work/g4eval/plots/crossing_t0%+.2f_%i_24022025.png", t0, runnumber));
      delete x;
      delete crossing;
    }

  //std::string filename = "collective_calibrations_24022025.csv";
  //std::string data = Form("%i,%i,%.8f,%.8f,%.8f,%f,%f,%f,%f,%i,%i",runnumber, nEntries,v_drift, new_drift, drift_err, t0_corr+t0, t0_err, misalignment, misalignment_err, nTotal, nZombies);
  //append_to_csv(filename, data);

  delete c;
  delete hist3D;
  delete fit;
  delete tree;
}

void drift_plotting_fit()
{

  drift_plotting_step(53285, 0, 0);
  return;

  std::map<int,float> drift_correction = csv_to_map("collective_calibrations_22022025.csv",2);
  std::map<int,float> t0_correction = csv_to_map("collective_calibrations_22022025.csv",5);
  for(const auto& [n, d] : drift_correction)
    {
      //drift_plotting_step(n, d, -5);
      drift_plotting_step(n, d, t0_correction[n]);
    }
}
