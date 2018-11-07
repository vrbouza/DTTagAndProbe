#include "DTTnPLocaltrigEff.h"

DTTnPLocaltrigEff::DTTnPLocaltrigEff(const std::string & configFile) : DTTnPBaseAnalysis(configFile)
{

}

void DTTnPLocaltrigEff::Loop()
{

  TFile outputFile(m_sampleConfig.outputFileName,"recreate");
  outputFile.cd();

  book();

  if (fChain == 0) return;

  Long64_t nentries = (m_sampleConfig.nEvents > 0 && 
		       fChain->GetEntriesFast() > m_sampleConfig.nEvents) ? 
    m_sampleConfig.nEvents : fChain->GetEntriesFast();

  std::cout << "[DTTnPLocaltrigEff::Loop] going to process "
	    << nentries << " entries\n" << std::flush;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry % 10000 == 0) 
	std::cout << "[DTTnPLocaltrigEff::Loop] processed : " 
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

void DTTnPLocaltrigEff::book()
{

  DTTnPBaseAnalysis::book();

  //  std::cout << " ** SM: booking trigger hostograms " << std::endl;

  m_plots["nOtherMatchedChVsEta"] = new TH2F("nOtherMatchedChVsEta",
					     "# of matched stations other than the one under investigation",
					     80,-1.2,1.2,5,-0.5,4.5);

  //  std::cout << " ** SM: DTTnPLocaltrigEff::book - booking trigger efficiency plots" << std::endl;

  for (Int_t iCh = 1; iCh < 5; ++iCh)
    {
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

      // ------------------ trigger efficiency (no segment requirement)

      hName = "trigeffAccVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance vs #eta;muon #eta;Efficiency",
				      80,-1.2,1.2);

    
      hName = "trigeffAccVsPhiPlus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance vs #phi for mu^{+};muon #phi;Efficiency",
				      80,-TMath::Pi(),TMath::Pi());

      hName = "trigeffAccVsPhiMinus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance vs #phi for mu^{-};muon #phi;Efficiency",
				      80,-TMath::Pi(),TMath::Pi());

      hName = "trigeffAccPhiVsEtaPlus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance #phi vs #eta for mu^{+};muon #phi;muon #eta",
				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAccPhiVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance #phi vs #eta;muon #phi;muon #eta",
				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);


      hName = "trigeffAccPhiVsEtaMinus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance #phi vs #eta for mu^{+};muon #phi;muon #eta",
				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAccVsPt" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance vs p_{T};muon p_{T};Efficiency",
				      50,0.,200.);

    

      hName = "trigeffAccVsLumi" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency x acceptance vs inst. lumi.; Efficiency",
				      25,0.,20000.);

      // ---

      hName = "trigeffVsPt" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency vs p_{T};muon p_{T};Efficiency",
				      50,0.,200.);

      hName = "trigeffVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency vs #eta;muon #eta;Efficiency",
				      80,-1.2,1.2);

      hName = "trigeffVsLumi" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency vs inst. lumi.;inst. lumi.;Efficiency",
				      25,0.,20000.);

      hName = "trigeffPhiVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
				      "trigger efficiency #phi vs #eta;muon #phi;muon #eta",
				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);
				     

      hName = "trigeffSecVsWh" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency sector vs wheel;sector;wheel",
      				      12,0.5,12.5,5,-2.5,2.5);
      hName = "trigeffChamb" + iChTag.str();
      m_plots[hName.c_str()] = new TH1F(hName.c_str(),
					"trigger efficiency chamber summary;efficiency;# chambers",
					400,0.,1.);

      // -------- efficiency when also a track segment in the station

     hName = "trigeffAcc_whensegm_VsRunNumb" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs RunNumber;Runnumber;Efficiency",
				      300,315000,318000);


      hName = "trigeffAcc_whensegm_VsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs #eta;muon #eta;Efficiency",
      				      80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_noZed_VsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm_noZed vs #eta;muon #eta;Efficiency",
      				      80,-1.2,1.2);
 
      hName = "trigeffAcc_whensegm_VsPhiPlus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs #phi for mu^{+};muon #phi;Efficiency",
      				      80,-TMath::Pi(),TMath::Pi());

      hName = "trigeffAcc_whensegm_VsPhiMinus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs #phi for mu^{-};muon #phi;Efficiency",
      				      80,-TMath::Pi(),TMath::Pi());

      hName = "trigeffAcc_whensegm_PhiVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm #phi vs #eta;muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_noZed_PhiVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm_noZed #phi vs #eta;muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_PhiVsEtaPlus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm #phi vs #eta for mu^{+};muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_PhiVsEtaPlus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm #phi vs #eta for mu^{+};muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_PhiVsEtaMinus" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm #phi vs #eta for mu^{+};muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);

      hName = "trigeffAcc_whensegm_VsPt" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs p_{T};muon p_{T};Efficiency",
      				      50,0.,200.);

      hName = "trigeffAcc_whensegm_SecVsWh" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm sector vs wheel;sector;wheel",
      				      12,0.5,12.5,5,-2.5,2.5);

      hName = "trigeffAcc_whensegm_noZed_SecVsWh" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm_noZed sector vs wheel;sector;wheel",
      				      12,0.5,12.5,5,-2.5,2.5);

      hName = "trigeffAcc_whensegm_Chamb" + iChTag.str();
      m_plots[hName.c_str()] = new TH1F(hName.c_str(),
      					"trigger efficiency acceptance whensegm chamber summary;efficiency;# chambers",
      					400,0.,1.);

      hName = "trigeffAcc_whensegm_VsLumi" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm vs inst. lumi.;inst. lumi.;Efficiency",
      				      25,0.,20000.);

      hName = "trigeffAcc_whensegm_noZed_VsLumi" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency x acceptance whensegm_noZed vs inst. lumi.;inst. lumi.;Efficiency",
      				      25,0.,20000.);
  
 
      hName = "trigeff_whensegm_VsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency whensegm vs #eta;muon #eta;Efficiency",
				      80,-1.2,1.2);

      hName = "trigeff_whensegm_VsPt" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency whensegm vs p_{T};muon p_{T};Efficiency",
      				      50,0.,200.);

      hName = "trigeff_whensegm_VsLumi" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency whensegm vs inst. lumi.;inst. lumi.;Efficiency",
      				      25,0.,20000.);
      
      hName = "trigeff_whensegm_PhiVsEta" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency whensegm #phi vs #eta;muon #phi;muon #eta",
      				      80,-TMath::Pi(),TMath::Pi(),80,-1.2,1.2);
				     

      hName = "trigeff_whensegm_SecVsWh" + iChTag.str();
      m_effs[hName] = new TEfficiency(hName.c_str(),
      				      "trigger efficiency whensegm sector vs wheel;sector;wheel",
      				      12,0.5,12.5,5,-2.5,2.5);

      hName = "trigeff_whensegm_Chamb" + iChTag.str();
      m_plots[hName.c_str()] = new TH1F(hName.c_str(),
      					"trigger efficiency whensegm chamber summary;efficiency;# chambers",
      					400,0.,1.);

      
      /// ---- Resolution plots 
//      hName = "trigRes_whensegm_SecVsWh" + iChTag.str();
//      m_plots[hName] = new TEfficiency(hName.c_str(),
//      				      "trigger Resolution whensegm sector vs wheel;sector;wheel",
//      				      12,0.5,12.5,5,-2.5,2.5);

      hName = "trigRes_whensegm_Chamb" + iChTag.str();
      m_plots[hName.c_str()] = new TH1F(hName.c_str(),
      					"trigger Resolution whensegm chamber summary;efficiency;# chambers",
      					400,0.,1.);

      hName = "trigRes_whensegm" + iChTag.str();
      m_plots[hName.c_str()] = new TH1F(hName.c_str(),
      					"trigger Resolution (phi) whensegm; 1 - l1 #phi / seg #phi;",
      					400,-2.,2.);
      
    }

 
}

void DTTnPLocaltrigEff::fill(const Int_t iMu)
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
         
      if ( nMatchInOtherCh >= m_tnpConfig.probe_minNMatchedSeg ) // ||
	//	   Mu_numberOfRPCLayers_rpc->at(iMu) >= m_tnpConfig.probe_minNRPCLayers )
	{
 
	  Int_t iPassingTMuxIn = getPassingTrig(iMu,iCh);  // TMuxOut primitive id in that station
          Int_t iPassingSeg = getPassingProbe(iMu,iCh); // track segment id in that station

          // ------  compute TMuxOut efficiency with no segment requirements

	  std::string hName = "trigeffAccVsEta" + iChTag.str();
	  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_eta->at(iMu));

	  hName = "trigeffAccVsLumi" + iChTag.str();
	  m_effs[hName]->Fill(iPassingTMuxIn >= 0,lumiperblock);

	  hName = "trigeffAccPhiVsEta" + iChTag.str();
	  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));

	  

		     

	  if (Mu_charge->at(iMu) == 1)
            {
              hName = "trigeffAccPhiVsEtaPlus" + iChTag.str();
              m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
            }

          if (Mu_charge->at(iMu) == -1)
            {
              hName = "trigeffAccPhiVsEtaMinus" + iChTag.str();
              m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
            }
	  // ------  compute TMuxIn efficiency with WITH segment requirements 
	 

	  bool segm_ok = false;
	  bool segm_ok_noZed = false;


	  segm_ok = (iPassingSeg>=0 && ((dtsegm4D_hasZed->at(iPassingSeg)>0 && iCh!=4 && dtsegm4D_phinhits->at(iPassingSeg)>=4) ||
					(iCh==4 && dtsegm4D_phinhits->at(iPassingSeg)>=4)));
	  segm_ok_noZed = (iPassingSeg>=0 && dtsegm4D_phinhits->at(iPassingSeg)>=4);
	  
	  // std::cout << " ** SM: segm_ok = " << segm_ok << std::endl;
	  
             
	  if(segm_ok_noZed) {
	    int sc=dtsegm4D_sector->at(iPassingSeg);
	    int wh=dtsegm4D_wheel->at(iPassingSeg);
	    
	    if(sc==13) {sc=4;}
	    else if(sc==14) {sc=10;}	
	    
	    std::string hName = "trigeffAcc_whensegm_noZed_VsEta" + iChTag.str();
	    m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_eta->at(iMu));
	    
	    hName = "trigeffAcc_whensegm_noZed_VsLumi" + iChTag.str();
	    m_effs[hName]->Fill(iPassingTMuxIn >= 0,lumiperblock);

	    hName = "trigeffAcc_whensegm_noZed_PhiVsEta" + iChTag.str();
	    m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
	    
	    hName = "trigeffAcc_whensegm_noZed_SecVsWh" + iChTag.str();
	    m_effs[hName]->Fill(iPassingTMuxIn >= 0,sc,wh);

	    if(segm_ok) {	

	      std::string hName = "trigeffAcc_whensegm_VsEta" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_eta->at(iMu));

	      // std::cout << " RUNnumber = " << runnumber << std::endl;
	      hName = "trigeffAcc_whensegm_VsRunNumb" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,runnumber);
	      
	      hName = "trigeffAcc_whensegm_VsLumi" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,lumiperblock);
	      
	      hName = "trigeffAcc_whensegm_PhiVsEta" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
	      
	      hName = "trigeffAcc_whensegm_SecVsWh" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,sc,wh);
	      

	      
	      if (Mu_charge->at(iMu) == 1)
		{
		  hName = "trigeffAcc_whensegm_PhiVsEtaPlus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
		}
	      
	      if (Mu_charge->at(iMu) == -1)
		{
		  hName = "trigeffAcc_whensegm_PhiVsEtaMinus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu),Mu_eta->at(iMu));
		}
	      
	    }
	  }
	  // -----------
	 
	  if (std::abs(probeVec.Eta()) < m_tnpConfig.probe_maxAbsEta[iCh - 1])
	    {

	      hName = "probePt" + iChTag.str(); 
	      m_plots[hName]->Fill(probeVec.Pt());

	      hName = "probeEta" + iChTag.str(); 
	      m_plots[hName]->Fill(probeVec.Eta());
	      
	      hName = "probePhi" + iChTag.str(); 
	      m_plots[hName]->Fill(probeVec.Phi());
	      
	      if (Mu_charge->at(iMu) == 1)
		{
		  hName = "trigeffAccVsPhiPlus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		}
	      
	      if (Mu_charge->at(iMu) == -1)
		{
		  hName = "trigeffAccVsPhiMinus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		}
	      
	      hName = "trigeffAccVsPt" + iChTag.str();
	      m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Pt());

	      if (Mu_charge->at(iMu) == 1)
		{
		  hName = "trigeffAccVsPhiPlus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		}
	      
	      if (Mu_charge->at(iMu) == -1)
		{
		  hName = "trigeffAccVsPhiMinus" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		}
	      
	     
	      // -----  compute TMuxIn efficiency with WITH segment requirements 
	     
              if(segm_ok_noZed) {

		int sc=dtsegm4D_sector->at(iPassingSeg);
		int wh=dtsegm4D_wheel->at(iPassingSeg);
		
		if(sc==13) {sc=4;}
		else if(sc==14) {sc=10;}
		  

		if(segm_ok) {

		

		  hName = "trigeffAcc_whensegm_VsPt" + iChTag.str();
		  m_effs[hName]->Fill(iPassingTMuxIn >= 0,probeVec.Pt());
		  
		  if (Mu_charge->at(iMu) == 1)
		    {
		      hName = "trigeffAcc_whensegm_VsPhiPlus" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		    }
		  
		  if (Mu_charge->at(iMu) == -1)
		    {
		      hName = "trigeffAcc_whensegm_VsPhiMinus" + iChTag.str();
		      m_effs[hName]->Fill(iPassingTMuxIn >= 0,Mu_phi->at(iMu));
		    }
		  
		}
	      }
       	      
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

  void DTTnPLocaltrigEff::endJob() 
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

  Int_t DTTnPLocaltrigEff::getPassingProbe(const Int_t iMu,
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

  Int_t DTTnPLocaltrigEff::getPassingTrig(const Int_t iMu,
					  const Int_t iCh) 
  {

    Int_t ibestTMuxIn = -1;
    Float_t bestDphi = 999.;	  
	  
    for (Int_t iMatch = 0; iMatch < Mu_nMatches->at(iMu); ++ iMatch)
      {
      

	Int_t whMu  = getXY<Int_t>(Mu_matches_Wh,iMu,iMatch);
	Int_t secMu = getXY<Int_t>(Mu_matches_Sec,iMu,iMatch);
	Int_t stMu  = getXY<Int_t>(Mu_matches_St,iMu,iMatch);
	// Float_t xMu = getXY<Float_t>(Mu_matches_x,iMu,iMatch); 
	// Float_t yMu = getXY<Float_t>(Mu_matches_y,iMu,iMatch); 
	Float_t phiMu = getXY<Float_t>(Mu_matches_phi,iMu,iMatch); 
       
      
	if (stMu == iCh)
	  {
	  
	    if(stMu==4 && secMu==13) {secMu=4;}
	    else if(stMu==4 && secMu==14 ) {secMu=10;}
	  
	    Int_t iTMuxIn = getPassingTrigInCh(iMu,stMu,secMu,
						whMu,phiMu);

	    if (iTMuxIn < 0)
	      continue;
	  
	    // Float_t xSeg = dtsegm4D_x_pos_loc->at(iSeg);
	    // Float_t ySeg = dtsegm4D_y_pos_loc->at(iSeg);

	    Float_t phiTrig = ltTwinMuxIn_phi->at(iTMuxIn);
	    phiTrig = phiTrig/4096 + ((3.1415927*30)/180)*(secMu-1);
	  
	  
	    // Float_t dX = std::abs(xSeg-xMu);
	    // Float_t dY = std::abs(ySeg-yMu);
	    Float_t Dphi = fabs(phiMu-phiTrig);
	    
	    if(Dphi < bestDphi)
	      {
		ibestTMuxIn = iTMuxIn;
		bestDphi = Dphi;
	      }
	    
	  }
      
      }
    if(ibestTMuxIn>=0) {
      
    }

    return ibestTMuxIn;

 

  }


  Int_t DTTnPLocaltrigEff::getPassingProbeInCh(const Int_t iMu,
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


  Int_t DTTnPLocaltrigEff::getPassingTrigInCh(const Int_t iMu,
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
    // if(ibestTMuxIn >=0 && bestTMuxInQual>=0) {
    //   std::cout << "----------------- " << std::endl;
    //   for (Int_t iTMuxIn = 0; iTMuxIn < NdtltTwinMuxIn; ++ iTMuxIn)
    //   	{
    // 	  std::cout << " **iTMuxIn ALL = " << iTMuxIn << std::endl;
    //   	  if(ltTwinMuxIn_wheel->at(ibestTMuxIn)==ltTwinMuxIn_wheel->at(iTMuxIn) &&
    //   	     ltTwinMuxIn_sector->at(ibestTMuxIn)==ltTwinMuxIn_sector->at(iTMuxIn) &&
    //   	     ltTwinMuxIn_station->at(ibestTMuxIn)==ltTwinMuxIn_station->at(iTMuxIn)
    // 	     )
	    
    // 	    // {std::cout << " iTMuxIn when wh st sec = " << iTMuxIn << "  is2nd = " <<  ltTwinMuxIn_is2nd->at(iTMuxIn) << std::endl;
	      
    // 	    //   if(ltTwinMuxIn_is2nd->at(iTMuxIn)==1) {
    // 	    // 	// compare phi, phib of the primitives
    // 	    // 	std::cout << " iTMuxIn = " << iTMuxIn << std::endl;
    // 	    // 	std::cout << " TMuxIn *** wh, st, sec = " << ltTwinMuxIn_wheel->at(ibestTMuxIn) << "," << ltTwinMuxIn_station->at(ibestTMuxIn) << "," << ltTwinMuxIn_sector->at(ibestTMuxIn) << " Rpc bit = " << ltTwinMuxIn_rpcbit->at(ibestTMuxIn) << std::endl;
    // 	    // 	std::cout << " phi best = " << ltTwinMuxIn_phi->at(ibestTMuxIn) << " phi is2nd = " << ltTwinMuxIn_phi->at(iTMuxIn) << std::endl;
    // 	    // 	std::cout << " phiB best = " << ltTwinMuxIn_phiB->at(ibestTMuxIn) << " phiB is2nd = " << ltTwinMuxIn_phiB->at(iTMuxIn) << std::endl; 
    // 	    // 	std::cout << " quality best = " << ltTwinMuxIn_quality->at(ibestTMuxIn) << "quality is2nd = " << ltTwinMuxIn_quality->at(iTMuxIn) << std::endl;
    // 	    //   }
    // 	    // }
	  
    // 	    }
      
    //   // for (Int_t ibmtfPh = 0; ibmtfPh < bmtfPhSize; ++ ibmtfPh)
    //   // 	{	  
    //   // 	  if(ltTwinMuxIn_wheel->at(ibestTMuxIn)==bmtfPhWh->at(ibmtfPh) &&
    //   // 	     ltTwinMuxIn_sector->at(ibestTMuxIn)==(bmtfPhSe->at(ibmtfPh)+1) &&
    //   // 	     ltTwinMuxIn_station->at(ibestTMuxIn)==bmtfPhSt->at(ibmtfPh) && 
    //   // 	     bmtfPhTs2Tag->at(ibmtfPh)==1)
    //   // 	    {
    //   // 		      // compare phi, phib of the primitives
    //   // 	      std::cout << " BMTF  wh, st, sec = " << ltTwinMuxIn_wheel->at(ibestTMuxIn) << "," << ltTwinMuxIn_station->at(ibestTMuxIn) << "," << ltTwinMuxIn_sector->at(ibestTMuxIn) + 1 << " is2nd = " <<  ltTwinMuxIn_rpcbit->at(ibestTMuxIn) << std::endl;
    //   // 	      std::cout << " phi best = " << ltTwinMuxIn_phi->at(ibestTMuxIn) << " bmtf phi is2nd = " << bmtfPhAng->at(ibmtfPh) << std::endl;
    //   // 	      std::cout << " phiB best = " << ltTwinMuxIn_phiB->at(ibestTMuxIn) << " bmtf phiB is2nd = " << bmtfPhBandAng->at(ibmtfPh) << std::endl; 
    //   // 	      std::cout << " quality best = " << ltTwinMuxIn_quality->at(ibestTMuxIn) << " bmtf quality is2nd = " << bmtfPhCode->at(ibmtfPh) << std::endl;      
    //   // 	    }
    //   // 	}
      
    // }
  
    
    return ibestTMuxIn;
    
    
  }
