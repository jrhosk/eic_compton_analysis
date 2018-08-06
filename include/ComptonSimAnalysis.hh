#ifndef compton_sim_analysis_hh
#define compton_sim_analysis_hh

#include <vector>
#include <string>

#include "TApplication.h"
#include "Generator.hh"
#include "vfTDCAnalysis.hh"

#include "boost/filesystem.hpp"

class ComptonSimAnalysis: public Generator, public std::ostream {

private:

public:


  static const double pressure_conversion = 1.315789e-12;
  static const int lvl_one_accept = 128;

  // Structure

  struct analysis {
    double luminosity;        // luminosity
    double cross_section;     // cross section
    double guassian_weight;   // weight to account for the fraction of the total guassian we are using
    double halo_ratio;        // halo ratio
    double halo_amplitude;
    double halo_scale_x;
    double halo_scale_y;
    double gaussian_weight;
    double aperture_size;
  };

  struct parameters { 

    double beam_energy;
    double laser_energy;
    double polarization;
    double sigma_ex;
    double sigma_ey;
    double sigma_g;
    double strip_number;
  };

  analysis compton;
  parameters beam;

  TApplication *app;

  std::string fFileLocation;
  std::string fFileRight;
  std::string fFileLeft;
  std::string fFileOutput;

  double fEnergyCut;       // Energy cut in MeV

  bool fGraphicsShow;
  bool fFileSet;
  bool fFileLeftSet;
  bool fFileRightSet;
  bool fComptonWeight;
  bool fBackgroundWeight;
  bool fHaloWeight;
  bool fApertureSize;
  bool fAsymmetryAnalysis;
  bool fVetrocAnalysis;
  bool fResolutionAnalysis;

  // Functions

  ComptonSimAnalysis();
  ~ComptonSimAnalysis();


  void GetSimulationFiles(const boost::filesystem::path& root, const std::string ext);
  void ReadConfiguration();
  void GetOptions(char **options);
  void ParseDetectorFile(std::fstream &file);
  void InitGraphicsEngine(int, char** );
  void RunGraphicsEngine();
  void AsymmetryAnalysis();
  void ScaleAsymmetry(TH1D *, TH1D *, TH1D *, int);

  double CalculateIntegratedCS();
  double CalculateIntegratedCS(double);
  double CalculateLuminosity(double);
  double CalculateGaussianWeight();
  double CalculateHaloFraction();

  bool OpenFile(TChain *);

};
#endif
