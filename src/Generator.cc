#define Generator_cxx
#include "../include/Generator.hh"
#include "TF2.h"
#include "TF1.h"
#include "TMath.h"
#include "../include/FontColor.hh"

#ifdef Generator_cxx

Generator::Generator(char *options): fGraphicsShow(false), fNumberEvents(1000)
{
  fFileName = (char *)"lund.dat";
  fPolarization = 1.0;
  sigma_x = 226.6e-4;
  sigma_y = 99e-4;

}

// Generator::Generator(double beam_e, double laser_e, double polar, char *options)
// {

//   beam_energy = beam_e;
//   laser_energy = laser_e;
//   polarization = polar;

// }

Generator::~Generator(){}


double Generator::CrossSection(double *x = 0, double *par = 0)
{

  double rho        = x[0];
  // double phi        = x[1];
  double b_energy   = par[0];
  double l_energy   = par[1];
  double P          = par[2];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  // Unpolarized cross section

  double term1 = (rho*rho*(1-alpha)*(1-alpha))/(1-rho*(1-alpha));
  double term2 = TMath::Power((1-rho*(1+alpha))/(1-rho*(1-alpha)), 2);
  double fdSig_dRho_0 = TMath::Power(electron_radius, 2)*alpha*(term1 + 1 + term2);

  // Polarized longitudinal cross section

  double term3 = (1-rho*(1+alpha))*(1-(1/TMath::Power( (1-rho*(1-alpha)),2) ));
  double fdSig_dRho_1 =TMath::Power(electron_radius, 2)*alpha*term3;

  // Total cross section for zero transverse polarization

  double fdSig_dRho = fdSig_dRho_0 + P*fdSig_dRho_1;

  return fdSig_dRho;

}

double Generator::BeamEnvelope(double *x, double *par){
  // Simple gaussian describing the core beam envelope.

  double pos = x[0];

  double sig = par[0];
  double diff_cross_section = 0;

  double top  = TMath::Power(pos, 2);
  double bot  = TMath::Power(sig, 2);

  diff_cross_section = TMath::Exp(-(top)/(2*bot) );


  return diff_cross_section;
}

void Generator::Initialize(){

  cs = new TF1("cs", this, &Generator::CrossSection, 0.0, 1.0, 3, "Generator", "CrossSection");
  cs->SetParameters(beam_energy, laser_energy, polarization);
  cs->SetNpx(1000);
  
}

void Generator::SetBeamEnergy(double energy){ 
  beam_energy = energy;
}
void Generator::SetLaserEnergy(double energy){ 
  laser_energy = energy;
}
void Generator::SetLaserWaveLength(double lambda){ 
  if(lambda > 0) laser_energy = h_planck*c_light/lambda;
}
void Generator::SetPolarization(double polar){ 
  polarization = polar;
}

double Generator::GetBeamEnergy(){return beam_energy;}

double Generator::GetLaserEnergy(){return laser_energy;}

double Generator::GetPolarization(){return polarization;}

double Generator::GetElectronMomentum(){return kinematics.electron_momentum;}

double Generator::GetPhotonMomentum(){return kinematics.photon_momentum;}

double Generator::GetElectronTheta(){return kinematics.electron_theta;}

double Generator::GetPhotonTheta(){return kinematics.photon_theta;}

double Generator::GetElectronPx(){return kinematics.px;}

double Generator::GetElectronPy(){return kinematics.py;}

double Generator::GetElectronPz(){return kinematics.pz;}

double Generator::GetPhotonPx(){return kinematics.kx;}

double Generator::GetPhotonPy(){return kinematics.ky;}

double Generator::GetPhotonPz(){return kinematics.kz;}


void Generator::CalculateKinematics()
{
  // double rho  = 0;
  double electron_prime = 0;

  kinematics.rho = cs->GetRandom();

  // kinematics.rho = rho;
  kinematics.photon_phi = gRandom->Uniform(0, 2*TMath::Pi()); // Sample from a uniform distribution of phi

  kinematics.kmax = 4*alpha*laser_energy*std::pow(beam_energy/electron_mass_c2, 2); // The maximum scattered photon energy or minimum electron energy

  kinematics.photon_momentum = kinematics.rho*kinematics.kmax;                               
  kinematics.electron_phi = -kinematics.photon_phi;                    

  electron_prime = beam_energy + laser_energy - kinematics.photon_momentum;

  // kinematics.electron_momentum = std::sqrt(std::pow(beam_energy, 2) - std::pow(electron_mass_c2, 2)); 
  kinematics.electron_momentum = std::sqrt(std::pow(electron_prime, 2) - std::pow(electron_mass_c2, 2)); 

  kinematics.photon_theta = std::sqrt(4*laser_energy*alpha/kinematics.photon_momentum-std::pow(electron_mass_c2/beam_energy, 2));
  // kinematics.photon_theta = std::sqrt( 4.*kinematics.photon_momentum/kinematics.kmax - 1./(alpha*std::pow(beam_energy/electron_mass_c2, 2)));
  kinematics.electron_theta = std::asin(kinematics.photon_momentum*std::sin(kinematics.photon_theta)/kinematics.electron_momentum); // check this

  kinematics.px = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::sin(kinematics.electron_phi);
  kinematics.py = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::cos(kinematics.electron_phi);
  kinematics.pz = kinematics.electron_momentum*std::cos(kinematics.electron_theta);

  kinematics.kx = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::sin(kinematics.photon_phi);
  kinematics.ky = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::cos(kinematics.photon_phi);
  kinematics.kz = kinematics.photon_momentum*std::cos(kinematics.photon_theta);

  kinematics.asymmetry = RhoToAsymmetry(beam_energy, laser_energy, kinematics.rho);
  // PrintAsymmetryInfo();
}

void Generator::PrintAsymmetryInfo()
{
  std::cout << kinematics.rho << "\t"
	    << kinematics.asymmetry
	    << std::endl;
}

void Generator::CalculateKinematics(event *kinematic){}

void Generator::GenerateAsymmetry(char* options)
{
  asym = new TF1("asym", this, &Generator::CalculateAsymmetry, 0.0, 1.0, 2, "Generator", "CalculateAsymmetry");
  asym->SetParameters(beam_energy, laser_energy);
  asym->SetNpx(1000);
}

double Generator::CalculateAsymmetry(double *x = 0, double *par = 0)
{

  double rho = x[0];
  double b_energy   = par[0];
  double l_energy   = par[1];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  double minus = rho*(1-alpha);
  double plus  = rho*(1+alpha);
  double term1 = 1/( (std::pow(minus, 2)/(1-minus)) + 1 + std::pow((1-plus)/(1-minus),2));
  double term2 = 1-plus;
  double term3 = 1-(1/std::pow(1-minus,2));

  double asymmetry = term1*term2*term3;

  return asymmetry;
}

double Generator::RhoToAsymmetry(double b_energy = 0, double l_energy = 0, double rho = 0)
{

  // double b_energy   = par[0];
  // double l_energy   = par[1];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  double minus = rho*(1-alpha);
  double plus  = rho*(1+alpha);
  double term1 = 1/( (std::pow(minus, 2)/(1-minus)) + 1 + std::pow((1-plus)/(1-minus),2));
  double term2 = 1-plus;
  double term3 = 1-(1/std::pow(1-minus,2));

  double asymmetry = term1*term2*term3;

  return asymmetry;
}

void Generator::OpenOutputFile()
{

  output.open(fFileName, std::fstream::out);

  if(!(output.is_open())){
    std::cerr << red << "Failure to open output file. Exiting." << white << std::endl;
    exit(1);
  }

}

void Generator::OpenOutputFile(char *filename)
{

  output.open(filename, std::fstream::out);

  if(!(output.is_open())){
    std::cerr << red << "Failure to open output file. Exiting." << white << std::endl;
    exit(1);
  }

}

void Generator::WriteHeader()
{
  event_counter++;

  output << "2 "
         << beam_energy << " "
         << laser_energy << " "
         << kinematics.kmax << " "
	 << polarization << " "
         // << "0. 0. 0. 0. 0. \n";
	 << kinematics.rho << " "
	 << kinematics.asymmetry << " "
	 << "0. 0. 0. \n";


}
void Generator::WriteEvent(int index, int pid, double px, double py, double pz, double momentum)
{

  output << index
         << " 0. 1 "
         << pid << " "
	 << " 0 0 "
	 // << kinematics.rho << " "
	 // << kinematics.asymmetry << " "
         << px << " "
         << py << " "
         << pz << " "
         << momentum << " "
    // << alpha << "-29.29464 0.0 -2287.855\n";                                                                                                                                 
	 << alpha << " "
         << kinematics.vx << " "
         << kinematics.vy << " "
         << kinematics.vz << "\n";
  
}

void Generator::ProcessEvent()
{

  WriteHeader(); // Write event header for a given event

  // Write event info for each final state particle.
  WriteEvent(1, pid_photon, -kinematics.kx, -kinematics.ky, -kinematics.kz, kinematics.photon_momentum);
  WriteEvent(2, pid_electron, -kinematics.px, -kinematics.py, -kinematics.pz, kinematics.electron_momentum);

  //  PrintEvent();

}

// void InitGeneratedAsymmetryGraph(int nevents)
// {

// }

void Generator::BuildGeneratedAsymmetryPlot()
{

}

void Generator::PrintEvent()
{

  std::cout << "\n=====================================\n" << std::endl;

  std::cout << "<<<< Electron \n" << std::endl;
  std::cout << "Electron theta: " << kinematics.electron_theta << std::endl;
  std::cout << "Electron phi: " << kinematics.electron_phi << std::endl;
  std::cout << "Electron Momentum: " << kinematics.electron_momentum << std::endl;
  std::cout << "     ---------------------------     " << std::endl;
  std::cout << "Electron px: " << kinematics.px << std::endl;
  std::cout << "Electron py: " << kinematics.py << std::endl;
  std::cout << "Electron pz: " << kinematics.pz << std::endl;

  std::cout << "\n\n\n" << std::endl;

  std::cout << "<<<< Photon \n" << std::endl;
  std::cout << "Photon theta: " << kinematics.photon_theta << std::endl;
  std::cout << "Photon phi: " << kinematics.photon_phi << std::endl;
  std::cout << "Photon Momentum: " << kinematics.photon_momentum << std::endl;
  std::cout << "     ---------------------------     " << std::endl;
  std::cout << "Photon px: " << kinematics.kx << std::endl;
  std::cout << "Photon py: " << kinematics.ky << std::endl;
  std::cout << "Photon pz: " << kinematics.kz << std::endl;
  std::cout << "\nPhoton max: " << kinematics.kmax << std::endl;

  std::cout << "\n=====================================\n" << std::endl;

}

void Generator::SetEventVertex(double x, double y, double z)
{

  kinematics.vx = x;    // Set scattering vertex coordinates (cm)
  kinematics.vy = y;
  kinematics.vz = z;

}

void Generator::CloseOutputFile()
{
  event_counter = 0;
  output.close();
}

int Generator::GetNumberEvents(){return fNumberEvents;}

TF1* Generator::GetFunction(char *option)
{

  std::string opt = std::string(option);

  if(strcmp(option, "cs") == 0)
    {
      return cs;
    }
    if(strcmp(option, "asym") == 0)
    {
      return asym;
    }
    return NULL;
}

void Generator::GetOptions(char **options)
{

  Int_t i = 0;

  TString flag;

  while(options[i] != NULL){
    flag = options[i];

    if(flag.CompareTo("--graphics", TString::kExact) == 0){
      flag.Clear();
      fGraphicsShow = true;

      std::cout << green << "<<<< Initializing TApplication for plots.\t" << white << std::endl;
    }
    if(flag.CompareTo("--filename", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fFileName = options[i + 1];
      std::cout << green << "<<<< Output file set to: " << fFileName << white << std::endl;
    }
    if(flag.CompareTo("--polarization", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fPolarization = atof(options[i + 1]);
      std::cout << green << "<<<< Polarization set to: " << fPolarization << white << std::endl;
    }
    if(flag.CompareTo("--energy", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      beam_energy = atof(options[i + 1]);
      std::cout << green << "<<<< Beam energy set to: " << beam_energy << white << std::endl;
    }
    if(flag.CompareTo("--sigmax", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      sigma_x = atof(options[i + 1]);
      std::cout << green << "<<<< Beam X-dispersion set to: " << sigma_x << white << std::endl;
    }
    if(flag.CompareTo("--sigmay", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      sigma_y = atof(options[i + 1]);
      std::cout << green << "<<<< Beam Y-dispersion set to: " << sigma_y << white << std::endl;
    }
    if(flag.CompareTo("--events", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fNumberEvents = atoi(options[i + 1]);
      std::cout << green << "<<<< Number of events generated: " << fNumberEvents << white << std::endl;
    }
    i++;
  }

}

void Generator::InitGraphicsEngine(int Argc, char **Argv)
{
  std::cout << green << "<<<< Initialize Graphics Engine." << white << std::endl;
  app = new TApplication("App", &Argc, Argv);

}

void Generator::RunGraphicsEngine()
{
  std::cout << green << "<<<< Running Graphics Engine." << white << std::endl;
  app->Run();
}

#endif
