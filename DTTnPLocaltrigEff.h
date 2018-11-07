#ifndef DTTnPLocaltrigEff_h
#define DTTnPLocaltrigEff_h

#include "DTTnPBaseAnalysis.h"

#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <bitset>
#include <regex>
#include <map>

class DTTnPLocaltrigEff : public DTTnPBaseAnalysis 
{

 public:
  DTTnPLocaltrigEff(const std::string & configFile);
  ~DTTnPLocaltrigEff() { };

  void Loop() override;
  
 protected:

  virtual void book() override;   
  virtual void fill(Int_t iMu) override;
  virtual void endJob() override;

  Int_t getPassingProbe(const Int_t iMu,
			const Int_t iCh);

  Int_t getPassingProbeInCh(const Int_t iMu,
                            const Int_t muSt,
                            const Int_t muSec,
                            const Int_t muWh,
                            const Int_t xMu,
                            const Int_t yMu);

  Int_t getPassingTrig(const Int_t iMu,
			const Int_t iCh);

  Int_t getPassingTrigInCh(const Int_t iMu,
                            const Int_t muSt,
                            const Int_t muSec,
                            const Int_t muWh,
                            const Float_t xMuphi);

};

#endif
