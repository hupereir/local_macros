void calc_scattering()
{

  double x0_lead = 0.5612; // cm
  double L_lead = 1; // cm

  double x0_iron = 1.757; // cn
  double L_iron = 1; // cm

  // double L_x0 = L_lead/x0_lead;
  // double L_x0 = L_iron/x0_iron;

  double L_x0 = 1.2/100;

  //
  double E = 4E3; // MeV

  // from https://iopscience.iop.org/article/10.1088/1742-6596/2349/1/012008
  double scattering = (13.6/E)*sqrt(L_x0)*(1.+0.0038*std::log(L_x0));
  double diffusion = scattering*0.4e6;

  std::cout << "scattering: " << scattering*1000 << " mrad" << std::endl;
  std::cout << "diffusion: " << diffusion << " um" << std::endl;




}
