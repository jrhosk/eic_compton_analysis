#ifndef vftdc_analysis_hh
#define vftdc_analysis_hh 

#include <vector>
#include <string>

#include "TApplication.h"
#include "ComptonSimAnalysis.hh"

#include "boost/filesystem.hpp"

class vfTDCAnalysis {

private:

  std::vector <int> channelList;

  static const int lvl_one_accept = 128;

public:

  bool fFileSet;

  std::string fFileLocation;

  vfTDCAnalysis();
  vfTDCAnalysis(std::string);
  ~vfTDCAnalysis();

  int vfTDCAnalyze();
  int GetChannelKey(int);

};
#endif
