#include "DTTnPLocalTrigRes.h"

// DT trigger
#include "DQM/DTMonitorModule/interface/DTTrigGeomUtils.h"

// Geometry
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"


DTTnPLocalTrigRes::DTTnPLocalTrigRes(const std::string & configFile) : DTTnPBaseAnalysis(configFile)
{

}

void DTTnPLocalTrigRes::Loop()
{

  TFile outputFile(m_sampleConfig.outputFileName,"recreate");
  outputFile.cd();

  book();

  if (fChain == 0) return;

  Long64_t nentries = (m_sampleConfig.nEvents > 0 && 
		       fChain->GetEntriesFast() > m_sampleConfig.nEvents) ? 
    m_sampleConfig.nEvents : fChain->GetEntriesFast();

  std::cout << "[DTTnPLocalTrigRes::Loop] going to process "
	    << nentries << " entries\n" << std::flush;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry % 10000 == 0) 
	std::cout << "[DTTnPLocalTrigRes::Loop] processed : " 
		  << jentry << " entries\r" << std::flush;

      bool hasGoodRun = false;

      for (const auto & run : m_sampleConfig.runs)
	{	  
	  if (run == 0 ||
	      run == runnumber)
	    {
	      hasGoodRun = true;
	      break;
	    }
	    
	}

      if(!hasGoodRun)
	continue;

      auto tnpPairs = tnpSelection();

      for(const auto & pair : tnpPairs) 
	{ 

	  fill(pair.second);

	}

    }

 
  
  endJob();

  outputFile.Write();
  outputFile.Close();

}

void DTTnPLocalTrigRes::book()
{

  DTTnPBaseAnalysis::book();

  
  m_plots["nOtherMatchedChVsEta"] = new TH2F("nOtherMatchedChVsEta",
					     "# of matched stations other than the one under investigation",
					     80,-1.2,1.2,5,-0.5,4.5);
  
  
  for (Int_t iCh = 1; iCh < 5; ++iCh) {
      std::stringstream iChTag;
      iChTag << "MB" << iCh;
      
      std::string hName = "probePt" + iChTag.str();
      m_plots[hName]  = new TH1F(hName.c_str(),
				 "probe p_{T};p_{T} [GeV];#entries/GeV",
				 200,0.,200.); 
      hName = "probeEta" + iChTag.str();
      m_plots[hName] = new TH1F(hName.c_str(),
				"probe #eta;#eta;#entries/0.05",
				56,-1.2,1.2); 
      hName = "probePhi" + iChTag.str();
      m_plots[hName] = new TH1F(hName.c_str(),
				"probe #phi;#phi;#entries/(pi*90)",
				180,-TMath::Pi(),TMath::Pi()); 
      
  }

  
  Int_t nPhiBins  = 401;
  Float_t rangePhi  = 10.025;
  Int_t nPhibBins = 401;
  Float_t rangePhiB = 10.025;
  
  /// ---- Resolution plots per sector, wheel & chamber: 
  for (Int_t iCh = 1; iCh <= 4; ++iCh) {
    std::stringstream station;  station << "_MB" << iCh;
    std::cout << station.str() << std::endl;
    
    std::string chTag = "_MB" + station.str();	  
    hName = "TM_PhiResidualIn_SecVsWh" + chTag.str();
    m_plots[hName.c_str()] = new TH2F(hName.c_str(), 
				      "Trigger local position In - Segment local position (correlated triggers); sector; wheel", 
      				      12,0.5,12.5,5,-2.5,2.5);
    
    hName = "TM_PhibResidualIn_SecVsWh" + chTag.str();
    m_plots[hName.c_str()] = new TH1F(hName.c_str(), ,
				      "Trigger local direction In - Segment local direction (correlated triggers); sector; wheel",
      				      12,0.5,12.5,5,-2.5,2.5);
    
    hName = "TM_PhiResidualOut_SecVsWh" + chTag.str();
    m_plots[hName.c_str()] = new TH1F(hName.c_str(),
				      "Trigger local position Out - Segment local position (correlated triggers); sector; wheel",
      				      12,0.5,12.5,5,-2.5,2.5);
    
    hName = "TM_PhibResidualOut_SecVsWh" + chTag.str();
    m_plots[hName.c_str()] = new TH1F(hName.c_str(),
				      "Trigger local direction Out - Segment local direction (correlated triggers); sector; wheel",
      				      12,0.5,12.5,5,-2.5,2.5);					

    
    for (Int_t iWh = -2; iWh <= 2; ++iWh)    {
      std::stringstream wheel; wheel << "_Wh" << iWh;
      std::cout << wheel.str() << std::endl;
      
      
	
      for (Int_t iSc = 1; iSc <= 12; ++iSc)    {
	std::stringstream sector; sector << "_Sec" << iWh;
	std::cout << sector.str() << std::endl;
	
	
	std::string chTag = "_W" + whell.str() + "_Sec" + sector.str() + "_MB" + station.str();	  
	hName = "TM_PhiResidualIn" + chTag.str();
	m_plots[hName.c_str()] = new TH1F(hName.c_str(), 
					  "Trigger local position In - Segment local position (correlated triggers); #phi", 
					  nPhiBins,-rangePhi,rangePhi);
	
	hName = "TM_PhibResidualIn" + chTag.str();
	m_plots[hName.c_str()] = new TH1F(hName.c_str(), ,
					  "Trigger local direction In - Segment local direction (correlated triggers); #phi",
					  nPhiBins,-rangePhi,rangePhi);
	
	hName = "TM_PhiResidualOut" + chTag.str();
	m_plots[hName.c_str()] = new TH1F(hName.c_str(),
					  "Trigger local position Out - Segment local position (correlated triggers); #phi",
					  nPhiBins,-rangePhi,rangePhi);
	
	hName = "TM_PhibResidualOut" + chTag.str();
	m_plots[hName.c_str()] = new TH1F(hName.c_str(),
					  "Trigger local direction Out - Segment local direction (correlated triggers); #phi",
					  nPhiBins,-rangePhi,rangePhi);
	
	hName = "TM_PhitkvsPhitrig" + chTag.str();
	m_plots[hName.c_str()] = new TH2F(hName.c_str(), 
					  "Local position: segment vs trigger",
					  100,-500.,500.,100,-500.,500.);
	
	hName = "TM_PhibtkvsPhibtrig" + chTag.str();
	m_plots[hName.c_str()] = new TH2F(hName.c_str(), 
					  "Local direction : segment vs trigger",
					  200,-40.,40.,200,-40.,40.);

	hName = "TM_PhibResidualvsTkPos" + chTag.str();
	m_plots[hName.c_str()] = new TH2F(hName.c_str(), 
					  "Local direction residual vs Segment Position",
					  100,-500.,500.,200,-10.,10.);
	
	hName = "TM_PhiResidualvsTkPos" + chTag.str();
	m_plots[hName.c_str()] = new TH2F(hName.c_str(), 
					  "Local Position residual vs Segment Position",
					  100,-500.,500.,200,-10.,10.);
      }
    }
  }
  
}

void DTTnPLocalTrigRes::fill(const Int_t iMu)
{

 

  TLorentzVector probeVec;
  probeVec.SetXYZM(Mu_px->at(iMu),
		   Mu_py->at(iMu),
		   Mu_pz->at(iMu),
		   0.106);

 

  for (Int_t iCh = 1; iCh < 5; ++iCh)
    {
      std::stringstream iChTag;
      iChTag << "MB" << iCh;
      
    
      Int_t nMatchInOtherCh = DTTnPBaseAnalysis::nMatchedCh(iMu,iCh);
     
      m_plots["nOtherMatchedChVsEta"]->Fill(Mu_eta->at(iMu),nMatchInOtherCh);
         
      if ( nMatchInOtherCh < m_tnpConfig.probe_minNMatchedSeg )  continue;
	  
      Int_t iPassingTMuxIn  = getPassingTrigIn(iMu,iCh);  // TMuxIn primitive id in that station
      Int_t iPassingTMuxOut = getPassingTrigOut(iMu,iCh);  // TMuxOut primitive id in that station
      Int_t iPassingSeg     = getPassingProbe(iMu,iCh); // track segment id in that station
      
		     
      // ------  compute TMuxIn efficiency with WITH segment requirements 
      bool segm_ok = false;
      segm_ok = (iPassingSeg>=0 && ((dtsegm4D_hasZed->at(iPassingSeg)>0 && iCh!=4 && dtsegm4D_phinhits->at(iPassingSeg)>=4) ||
				    (iCh==4 && dtsegm4D_phinhits->at(iPassingSeg)>=4)));
      
      if (!segm_ok) continue;
      
      // -----------
      if (std::abs(probeVec.Eta()) > m_tnpConfig.probe_maxAbsEta[iCh - 1]) continue;
      
      for (Int_t iMatch = 0; iMatch < Mu_nMatches->at(iMu); ++iMatch)   // this is in full acceptance region (out of borders)
	{
	  Int_t whMu  = getXY<Int_t>(Mu_matches_Wh,iMu,iMatch);
	  Int_t secMu = getXY<Int_t>(Mu_matches_Sec,iMu,iMatch);
	  Int_t stMu  = getXY<Int_t>(Mu_matches_St,iMu,iMatch);
	  Float_t xMu = getXY<Float_t>(Mu_matches_x,iMu,iMatch);
	  Float_t yMu = getXY<Float_t>(Mu_matches_y,iMu,iMatch);
	  Float_t phiMu = getXY<Float_t>(Mu_matches_phi,iMu,iMatch);
	  
	  Float_t xBorderMu = getXY<Float_t>(Mu_matches_edgeX,iMu,iMatch);
	  Float_t yBorderMu = getXY<Float_t>(Mu_matches_edgeY,iMu,iMatch);
	  
	  if (stMu == iCh &&
	      xBorderMu < m_tnpConfig.probe_maxBorderDx &&
	      yBorderMu < m_tnpConfig.probe_maxBorderDy )
		    { 
		      
		      int tr_secMu;
		      tr_secMu=secMu;
		      if(stMu==4 && secMu==13) {tr_secMu=4; }
		      else if(stMu==4 && secMu==14) {tr_secMu=10; }

		      iPassingTMuxIn =  getPassingTrigInCh(iMu,stMu,tr_secMu,whMu,phiMu);
		      iPassingSeg =  getPassingProbeInCh(iMu,stMu,secMu,whMu,xMu,yMu);

		      segm_ok = false;
		      segm_ok = (iPassingSeg>=0 && ((dtsegm4D_hasZed->at(iPassingSeg)>0 && iCh!=4 && dtsegm4D_phinhits->at(iPassingSeg)>=4) ||
						    (iCh==4 && dtsegm4D_phinhits->at(iPassingSeg)>=4)));
		      segm_ok_noZed = false;
		      segm_ok_noZed = (iPassingSeg>=0 && dtsegm4D_phinhits->at(iPassingSeg)>=4);
		      
		      
		      
		      hName = "trigeffVsPt" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Pt());
		      
		      std::string hName = "trigeffVsEta" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Eta());
		      
		      hName = "trigeffPhiVsEta" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Phi(),probeVec.Eta());
		      
		      hName = "trigeffSecVsWh" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,tr_secMu,whMu);
		      
		      
		      
		      hName = "trigeffVsLumi" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,lumiperblock);
		      
		      
		      if(segm_ok_noZed) {

			if(segm_ok) {
			  hName = "trigeff_whensegm_VsPt" + iChTag.str();
			  m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Pt());
			  
			  hName = "trigeff_whensegm_PhiVsEta" + iChTag.str();
			  m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Phi(),probeVec.Eta());
			  
			  hName = "trigeff_whensegm_VsEta" + iChTag.str();
			  m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Eta());
			  
			  
			  hName = "trigeff_whensegm_SecVsWh" + iChTag.str();
			  m_effs[hName]->Fill(iPassingTMuxIn >= 0,tr_secMu,whMu);
			  
			  hName = "trigeff_whensegm_VsLumi" + iChTag.str();
			  m_effs[hName]->Fill(iPassingTMuxIn >= 0,lumiperblock);
			  
			  
			  hName = "trigRes_whensegm" + iChTag.str();
			  if (iPassingSeg >=0 && iPassingTMuxIn>=0) {
			    Float_t muphi = dtsegm4D_phi->at(iPassingSeg);
			    Float_t trigphi = ltTwinMuxIn_phi->at(iPassingTMuxIn);
			    trigphi = trigphi/4096 + ((3.1415927*30)/180)*(secMu-1);
			    m_plots[hName]->Fill(1.-trigphi/muphi);
			  }
			}
		      }
		    }
		}
	    }
	}
    }   
}

  void DTTnPLocalTrigRes::endJob() 
  {
    Int_t iS;
    for (Int_t iCh = 1; iCh < 5; ++iCh)
      {
	std::stringstream iChTag;
	iChTag << "MB" << iCh;

	for (Int_t iWh = 1; iWh < 6; ++iWh)
	  {

	    for (Int_t iSec = 1; iSec < 15; ++iSec)
	      {

		if (iSec > 12 && iCh !=4)
		  continue;
		iS=iSec;
		if(iSec==13 && iCh==4) {iS=4;}
		else if(iSec==14 && iCh==4) {iS=10;}
	      
		std::string hName   = "trigeffChamb" + iChTag.str();
		std::string effName = "trigeffSecVsWh" + iChTag.str();
		m_plots[hName.c_str()]->Fill(m_effs[effName.c_str()]->GetEfficiency(m_effs[effName.c_str()]->GetGlobalBin(iS,iWh)));

		hName   = "trigeff_whensegm_Chamb" + iChTag.str();
		effName = "trigeff_whensegm_SecVsWh" + iChTag.str();
		m_plots[hName.c_str()]->Fill(m_effs[effName.c_str()]->GetEfficiency(m_effs[effName.c_str()]->GetGlobalBin(iS,iWh)));

		hName   = "trigeffAcc_whensegm_Chamb" + iChTag.str();
		effName = "trigeffAcc_whensegm_SecVsWh" + iChTag.str();
		m_plots[hName.c_str()]->Fill(m_effs[effName.c_str()]->GetEfficiency(m_effs[effName.c_str()]->GetGlobalBin(iS,iWh)));
		
//		hName = "trigRes_whensegm_Cham" + iChTag.str();
//		HistoName = "trigRes_whensegm_SecVsWh"  + iChTag.str();
//		m_plots[hName.c_str()]->Fill(m_plots[HistoName.c_str()]->GetBinContent(m_plots[HistoName.c_str()]->GetBin(iS,iWh));
		// hName   = "trigeffAcc_whensegm_noZed_Chamb" + iChTag.str();
		// effName = "trigeffAcc_whensegm_NoZed_SecVsWh" + iChTag.str();
		// m_plots[hName.c_str()]->Fill(m_effs[effName.c_str()]->GetEfficiency(m_effs[effName.c_str()]->GetGlobalBin(iS,iWh)));
	      
	      }
	  }
      }
  }

  Int_t DTTnPLocalTrigRes::getPassingProbe(const Int_t iMu,
					   const Int_t iCh) 
  {

    Int_t iBestSeg   = -1;
    Float_t bestSegDr = 999.;	  
	  
    for (Int_t iMatch = 0; iMatch < Mu_nMatches->at(iMu); ++ iMatch)
      {
	Int_t whMu  = getXY<Int_t>(Mu_matches_Wh,iMu,iMatch);
	Int_t secMu = getXY<Int_t>(Mu_matches_Sec,iMu,iMatch);
	Int_t stMu  = getXY<Int_t>(Mu_matches_St,iMu,iMatch);
	Float_t xMu = getXY<Float_t>(Mu_matches_x,iMu,iMatch); 
	Float_t yMu = getXY<Float_t>(Mu_matches_y,iMu,iMatch); 
      

	if (stMu == iCh)
	  {
	  
	    Int_t iSeg = getPassingProbeInCh(iMu,stMu,secMu,
					     whMu,xMu,yMu);

	    if (iSeg < 0)
	      continue;
	  
	    Float_t xSeg = dtsegm4D_x_pos_loc->at(iSeg);
	    Float_t ySeg = dtsegm4D_y_pos_loc->at(iSeg);
	  
	    Float_t dX = std::abs(xSeg-xMu);
	    Float_t dY = std::abs(ySeg-yMu);
	    Float_t dR = sqrt(dX*dX + dY*dY);
	  
	    if(dR < bestSegDr)
	      {
		iBestSeg = iSeg;
		bestSegDr = dR;
	      }
	  
	  }
      
      }

    return iBestSeg;

  }

Int_t DTTnPLocalTrigRes::getPassingTrigIn(const Int_t iMu,
					  const Int_t iCh) 
{
  
  Int_t ibestTMuxIn = -1;
  Float_t bestDphi = 999.;	  
  
  for (Int_t iMatch = 0; iMatch < Mu_nMatches->at(iMu); ++ iMatch)      {
    
    Int_t whMu  = getXY<Int_t>(Mu_matches_Wh,iMu,iMatch);
    Int_t secMu = getXY<Int_t>(Mu_matches_Sec,iMu,iMatch);
    Int_t stMu  = getXY<Int_t>(Mu_matches_St,iMu,iMatch);
    Float_t phiMu = getXY<Float_t>(Mu_matches_phi,iMu,iMatch); 
       
    if (stMu == iCh)
      {
	if(stMu==4 && secMu==13) {secMu=4;}
	else if(stMu==4 && secMu==14 ) {secMu=10;}
	
	Int_t iTMuxIn = getPassingTrigInCh(iMu,stMu,secMu,whMu,phiMu);
	if (iTMuxIn < 0)
	  continue;
	

	Float_t phiTrig = ltTwinMuxIn_phi->at(iTMuxIn);
	phiTrig = phiTrig/4096 + ((3.1415927*30)/180)*(secMu-1);
	  
	Float_t Dphi = fabs(phiMu-phiTrig);
	if(Dphi < bestDphi)
	  {
	    ibestTMuxIn = iTMuxIn;
	    bestDphi = Dphi;
	  }
	
      }
    
  }
  
  return ibestTMuxIn;

}

Int_t DTTnPLocalTrigRes::getPassingTrigOut(const Int_t iMu,
					   const Int_t iCh) 
  {
    
    Int_t ibestTMuxOut = -1;
    Float_t bestDphi = 999.;	  
    
    for (Int_t iMatch = 0; iMatch < Mu_nMatches->at(iMu); ++ iMatch)
      {

	Int_t whMu  = getXY<Int_t>(Mu_matches_Wh,iMu,iMatch);
	Int_t secMu = getXY<Int_t>(Mu_matches_Sec,iMu,iMatch);
	Int_t stMu  = getXY<Int_t>(Mu_matches_St,iMu,iMatch);
	Float_t phiMu = getXY<Float_t>(Mu_matches_phi,iMu,iMatch); 
       
      	if (stMu == iCh)
	  {
	    
	    if(stMu==4 && secMu==13) {secMu=4;}
	    else if(stMu==4 && secMu==14 ) {secMu=10;}
	  
	    Int_t iTMuxOut = getPassingTrigOutCh(iMu,stMu,secMu,
						 whMu,phiMu);

	    if (iTMuxOut < 0)
	      continue;
	  
	    Float_t phiTrig = ltTwinMuxOut_phi->at(iTMuxIn);
	    phiTrig = phiTrig/4096 + ((3.1415927*30)/180)*(secMu-1);
	  
	    Float_t Dphi = fabs(phiMu-phiTrig);
	    
	    if(Dphi < bestDphi)
	      {
		ibestTMuxOut = iTMuxOut;
		bestDphi = Dphi;
	      }
	    
	  }
      
      }
    return ibestTMuxOut;
  }


Int_t DTTnPLocalTrigRes::getPassingProbeInCh(const Int_t iMu,
					     const Int_t stMu,
					     const Int_t secMu,
					     const Int_t whMu,
					     const Int_t xMu,
					     const Int_t yMu)
{
  
    Int_t iBestSeg   = -1;
    Float_t bestSegDr = 999.;
  
    for (Int_t iSeg = 0; iSeg < Ndtsegments; ++ iSeg)
      {
      
	Int_t whSeg  = dtsegm4D_wheel->at(iSeg);
	Int_t secSeg = dtsegm4D_sector->at(iSeg);
	Int_t stSeg  = dtsegm4D_station->at(iSeg);
	Float_t xSeg = dtsegm4D_x_pos_loc->at(iSeg);
	Float_t ySeg = dtsegm4D_y_pos_loc->at(iSeg);
      
	Float_t dX = std::abs(xSeg-xMu);
	Float_t dY = std::abs(ySeg-yMu);
	Float_t dR = sqrt(dX*dX + dY*dY);

	// if(stSeg ==4 && secSeg==13)  {secSeg=4;}
	// else if(stSeg ==4 && secSeg==14)  {secSeg=10;}
      
	if(whMu  == whSeg  &&
	   secMu == secSeg &&
	   stMu  == stSeg  &&
	   dR < bestSegDr  &&
	   dR < m_tnpConfig.passing_probe_maxTkSegDr &&
	   dX < m_tnpConfig.passing_probe_maxTkSegDx &&
	   dY < m_tnpConfig.passing_probe_maxTkSegDy )
	  {
	    iBestSeg = iSeg;
	    bestSegDr = dR;
	  }
      }

    return iBestSeg;

  }


  Int_t DTTnPLocalTrigRes::getPassingTrigInCh(const Int_t iMu,
					      const Int_t stMu,
					      const Int_t secMu,
					      const Int_t whMu,
					      const Float_t phiMu
					      )
  {
  
    Int_t ibestTMuxIn   = -1;
    Float_t bestTMuxInQual = -1.;
    //       std::cout << " ******* another call to local trig eff " << std::endl;
    for (Int_t iTMuxIn = 0; iTMuxIn < NdtltTwinMuxIn; ++ iTMuxIn)
      {
      
	Int_t whTMuxIn  = ltTwinMuxIn_wheel->at(iTMuxIn);
	Int_t secTMuxIn = ltTwinMuxIn_sector->at(iTMuxIn);
	Int_t stTMuxIn  = ltTwinMuxIn_station->at(iTMuxIn);
	Int_t bxTMuxIn  = ltTwinMuxIn_bx->at(iTMuxIn);
	Int_t phiTMuxIn  = ltTwinMuxIn_phi->at(iTMuxIn);
	Int_t phiBTMuxIn  = ltTwinMuxIn_phiB->at(iTMuxIn);
	Int_t qualityTMuxIn  = ltTwinMuxIn_quality->at(iTMuxIn);

	if(stTMuxIn ==4 && secTMuxIn==13)  {secTMuxIn=4;}
	else if(stTMuxIn ==4 && secTMuxIn==14)  {secTMuxIn=10;}
      
	if(whMu  == whTMuxIn   &&
	   secMu == secTMuxIn  &&
	   stMu  == stTMuxIn   && 
	   bxTMuxIn == 0       &&
	   qualityTMuxIn > bestTMuxInQual
	   )
	
	  {
	    ibestTMuxIn = iTMuxIn;
	    bestTMuxInQual = qualityTMuxIn;	  
	  }
      }
  
    
    return ibestTMuxIn;
    
    
  }


  Int_t DTTnPLocalTrigRes::getPassingTrigOutCh(const Int_t iMu,
					      const Int_t stMu,
					      const Int_t secMu,
					      const Int_t whMu,
					      const Float_t phiMu
					      )
  {
  
    Int_t ibestTMuxOut   = -1;
    Float_t bestTMuxOutQual = -1.;
    //       std::cout << " ******* another call to local trig eff " << std::endl;
    for (Int_t iTMuxOut = 0; iTMuxOut < NdtltTwinMuxOut; ++ iTMuxOut)
      {
      
	Int_t whTMuxOut  = ltTwinMuxOut_wheel->at(iTMuxOut);
	Int_t secTMuxOut = ltTwinMuxOut_sector->at(iTMuxOut);
	Int_t stTMuxOut  = ltTwinMuxOut_station->at(iTMuxOut);
	Int_t bxTMuxOut  = ltTwinMuxOut_bx->at(iTMuxOut);
	Int_t phiTMuxOut  = ltTwinMuxOut_phi->at(iTMuxOut);
	Int_t phiBTMuxOut  = ltTwinMuxOut_phiB->at(iTMuxOut);
	Int_t qualityTMuxOut  = ltTwinMuxOut_quality->at(iTMuxOut);

	if(stTMuxOut ==4 && secTMuxOut==13)  {secTMuxOut=4;}
	else if(stTMuxOut ==4 && secTMuxOut==14)  {secTMuxOut=10;}
      
	if(whMu  == whTMuxOut   &&
	   secMu == secTMuxOut  &&
	   stMu  == stTMuxOut   && 
	   bxTMuxOut == 0       &&
	   qualityTMuxOut > bestTMuxOutQual
	   )
	
	  {
	    ibestTMuxOut = iTMuxOut;
	    bestTMuxOutQual = qualityTMuxOut;	  
	  }
      }
  
    
    return ibestTMuxOut;
    
    
  }

void DTTnPLocalTrigRes::computeSCCoordinates(float localDirX, 
					     float localDirY, 
					     float localDirZ, 
					     int sector, int station,
					     int& scsec, float& x, float& xdir, float& y, float& ydir)
{
  xdir = TMath::ATan(localDirX / localDirY)*TMath::RadToDeg();
  ydir = TMath::ATan(localDirY / localDirZ)*TMath::RadToDeg();
  
  scsec = sector>12 ? sector==13 ? 4 : 10 : sector;
  float xcenter = (scsec==4||scsec==10) ? (sector-12.9)/abs(sector-12.9)*xCenter_[(sector==10||sector==14)] : 0.;
  x = localPosX+xcenter*(station==4);
  y = localPosY;

}
