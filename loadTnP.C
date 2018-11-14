#include <string>
#include "TROOT.h"

void loadSegmentTnP()
{

  gROOT->ProcessLine(".L DTAnalyzer.C++");
  gROOT->ProcessLine(".L DTTnPConfig.C++");
  gROOT->ProcessLine(".L DTTnPBaseAnalysis.C++");
  gROOT->ProcessLine(".L DTTnPSegmentEff.C++");

}

void loadTMuxOutTnP()
{

  gROOT->ProcessLine(".L DTAnalyzer.C++");
  gROOT->ProcessLine(".L DTTnPConfig.C++");
  gROOT->ProcessLine(".L DTTnPBaseAnalysis.C++");
  gROOT->ProcessLine(".L DTTnPLocaltrigEff.C++");

}
void loadTPTnP()
{

  gROOT->ProcessLine(".L DTAnalyzer.C++");
  gROOT->ProcessLine(".L DTTnPConfig.C++");
  gROOT->ProcessLine(".L DTTnPBaseAnalysis.C++");
  gROOT->ProcessLine(".L DTTnPLocaltrigEff.C++");
  gROOT->ProcessLine(".L DTTnPLocalTrigRes.C++");
  
}

