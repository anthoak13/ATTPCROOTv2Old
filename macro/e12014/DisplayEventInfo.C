void DisplayEventInfo()
{
	return 0;
}
	
void displayFilteredData(int tpcRun = 206)
{

   TChain tpc_tree("cbmsim");
   tpc_tree.Add(TString::Format("/mnt/analysis/e12014/TPC/filterTesting/run_%04d.root", tpcRun));

   TTreeReader reader(&tpc_tree);
   TTreeReaderValue<TClonesArray> event(reader, "AtRawEvent");
   TTreeReaderValue<TClonesArray> eventFilter(reader, "AtRawEventFiltered");

   const int numEvents = 2;

   
   // Loop through every event
   for (int i = 0; i < numEvents; ++i) {
      if (!reader.Next())
         break;

      // Get the unfiltered and filtered events
      AtRawEvent *eventPtr = (AtRawEvent *)(event->At(0));
      AtRawEvent *filterPtr = (AtRawEvent *)(eventFilter->At(0));	 
	  
      // Get the pad to display on both
      auto pad = eventPtr->GetPad(0);
      auto padNum = pad->GetPadNum();

      // Pointer to hold the filtered pad to display
      AtPad *filterPad;

      // Loop through and look for filtered pad.
      for (int i = 0; i < filterPtr->GetNumPads(); ++i) {
         filterPad = filterPtr->GetPad(i);
         if (filterPad->GetPadNum() == padNum)
            break;
      } // end loop over pads

      // Show the two pads
      TH1F *raw = new TH1F("rawHist", "RawData", 512, 0, 511);
      TH1F *filtered = new TH1F("filteredHist", "FilteredData", 512, 0, 511);
	  double maxADC = 0;
      for (int i = 0; i < 512; ++i) {
         raw->SetBinContent(i, pad->GetADC(i));
         filtered->SetBinContent(i, filterPad->GetADC(i));
		 
		 if (pad->GetADC(i) > maxADC){
			 maxADC = pad->GetADC(i);
		 }
		 
      }
      TCanvas *c1 = new TCanvas();
      raw->Draw();
      TCanvas *c3 = new TCanvas();
      filtered->Draw();
	  
	  int numberOfPads = eventPtr->GetNumPads();

	  
		std::cout << TString::Format("For Event %d" ,i) << std::endl;
		std::cout << TString::Format("Total pads hit = %d, Max value in shown pad = %f", numberOfPads, maxADC)<< std::endl; 
		std::cout << "\n" << std::endl;
   }
}

void displayEventData(int tpcRun = 206)
{

   TChain tpc_tree("cbmsim");
   tpc_tree.Add(TString::Format("/mnt/analysis/hira_collaboration/e12014/Joe/UnpackedRuns/run_%04d.root", tpcRun));

   TTreeReader reader(&tpc_tree);
   TTreeReaderValue<TClonesArray> event(reader, "AtRawEvent");
   TTreeReaderValue<TClonesArray> eventH(reader, "AtEventH");

   const int numEvents = 200;
   int FissionEvents = 0;

   // Loop through every event
   for (int i = 0; i < numEvents; ++i) {
      if (!reader.Next())
         break;

	double maxADC = 0; 
	 
      // Get the events
      AtRawEvent *eventPtr = (AtRawEvent *)(event->At(0));
	//  AtEventH *eventHPtr = (AtEventH *)(eventH->At(0));		 

	
	  // Get the pad to display on both
      auto pad = eventPtr->GetPad(0);
      auto padNum = pad->GetPadNum();

      // Pointer to hold the filtered pad to display
      AtPad *filterPad;	  
	  
      // Loop through pads in event.
      for (int i = 0; i < eventPtr->GetNumPads(); ++i) {
        pad = eventPtr->GetPad(i);
		for (int j = 0; j < 512; ++j) {
		 if (pad->GetADC(j) > maxADC){
			 maxADC = pad->GetADC(j);
			 filterPad = eventPtr->GetPad(i);
		 } 
		}
		
      } // end loop over pads

      // Show the two pads
//      TH1F *raw = new TH1F(TString::Format("rawHist%d", (i+1)), TString::Format("rawData%d", (i+1)), 512, 0, 511);
//      for (int i = 0; i < 512; ++i) { 
//         raw->SetBinContent(i, filterPad->GetADC(i));
//      }
//      TCanvas *c1 = new TCanvas();
//      raw->Draw();
	  
	  	  
	  int numberOfPads = eventPtr->GetNumPads();
	  if(numberOfPads>300){FissionEvents++;}
//	  double totalCharge = eventHPtr->GetQHit();
	  
		std::cout << TString::Format("For Event %d" ,i) << std::endl;
		std::cout << TString::Format("Total pads hit = %d, Max ADC in any pad = %f", numberOfPads, maxADC)<< std::endl; 
		std::cout << "" << std::endl;
   }
      std::cout << TString::Format("Possible Fission Events (by Pad Hits) = %d" , FissionEvents) << std::endl;
}

void displayTotalCharge(int tpcRun = 206, int tpcEvent = 80)
{

   TChain tpc_tree("cbmsim");
   tpc_tree.Add(TString::Format("/mnt/analysis/hira_collaboration/e12014/Joe/UnpackedRuns/run_%04d.root", tpcRun));

   TTreeReader reader(&tpc_tree);
   TTreeReaderValue<TClonesArray> event(reader, "AtRawEvent");
   TTreeReaderValue<TClonesArray> eventH(reader, "AtEventH");

   reader.SetEntry(tpcEvent);
   
	double maxADC = 0; 
	 
      // Get the events
      AtRawEvent *eventPtr = (AtRawEvent *)(event->At(0));
	  AtEvent *eventHPtr = (AtEvent *)(eventH->At(0));		 

	
	  // Get the pad to display on both
      auto pad = eventPtr->GetPad(0);
      auto padNum = pad->GetPadNum();

      // Pointer to hold the filtered pad to display
      AtPad *filterPad;	  
	  
      // Loop through pads in event.
      for (int i = 0; i < eventPtr->GetNumPads(); ++i) {
        pad = eventPtr->GetPad(i);
		for (int j = 0; j < 512; ++j) {
		 if (pad->GetADC(j) > maxADC){
			 maxADC = pad->GetADC(j);
			 filterPad = eventPtr->GetPad(i);
		 } 
		}
		
      } // end loop over pads

      // Show the two pads
      TH1F *adc = new TH1F(TString::Format("adcHist%d", (tpcEvent)), TString::Format("adcData%d", (tpcEvent)), 512, 0, 511);
	  TH1F *charge = new TH1F(TString::Format("chargeHist%d", (tpcEvent)), TString::Format("chargeData%d", (tpcEvent)), 512, 0, 511);
      for (int i = 0; i < 512; ++i) { 
         adc->SetBinContent(i, filterPad->GetADC(i));
      }
      TCanvas *c1 = new TCanvas();
      adc->Draw();
	  
	  	  
	  int numberOfPads = eventPtr->GetNumPads();
	  double totalCharge = eventHPtr->GetEventCharge();
	  
		std::cout << TString::Format("For Event %d" , tpcEvent) << std::endl;
		std::cout << TString::Format("Total pads hit = %d, Max ADC in any pad = %f, Total Charge = %f", numberOfPads, maxADC, totalCharge)<< std::endl; 
		std::cout << "" << std::endl;
}

void plotInfoVsEvent(int tpcRun = 206)
{
	
   TChain tpc_tree("cbmsim");
   tpc_tree.Add(TString::Format("/mnt/analysis/hira_collaboration/e12014/Joe/UnpackedRuns/run_%04d.root", tpcRun));

   TTreeReader reader(&tpc_tree);
   TTreeReaderValue<TClonesArray> event(reader, "AtRawEvent");
   TTreeReaderValue<TClonesArray> eventH(reader, "AtEventH");

	 


   const int numEvents = 1000;
   int FissionEvents = 0; //To count events that pass certain criteria
   int indexEvents[200]; //To plot against event number
   double maxADCs[200]; //To store max ADC signal in any pad for each event
   double totalCharges[200]; //To store total charge of each event
   int totalPadsHit[200]; //To store number of hit pads for each event

   // Loop through every event
   for (int i = 0; i < numEvents; ++i) {
      if (!reader.Next())
         break;
	 
	 
   // Get the events
   AtRawEvent *eventPtr = (AtRawEvent *)(event->At(0));
   AtEvent *eventHPtr = (AtEvent *)(eventH->At(0));	
	double maxADC = 0; 

	  // Get the pad to display on both
      auto pad = eventPtr->GetPad(0);
	  
      // Loop through pads in event.
      for (int k = 0; k < eventPtr->GetNumPads(); ++k) 
	  {
        pad = eventPtr->GetPad(k);
		for (int j = 0; j < 512; ++j) 
		{
		 if (pad->GetADC(j) > maxADC)
		 {
			 maxADC = pad->GetADC(j);
		 } 
		}
      } // end loop over pads

	  	  
		indexEvents[i] = i;
		totalCharges[i] = eventHPtr->GetEventCharge();
		totalPadsHit[i] = eventPtr->GetNumPads();
		maxADCs[i] = maxADC;
		int numberOfPads = eventPtr->GetNumPads();
		if(numberOfPads>300){FissionEvents++;}
   }
   
   
      std::cout << TString::Format("Possible Fission Events (by Pad Hits) = %d" , FissionEvents) << std::endl;
	  
	  
	  
	  TH1F *adc = new TH1F("adcHist",TString::Format("Max ADC Value (Any Pad) Vs. Event Number, Run %d", tpcRun), 200, 0, 199);
	  TH1F *numpads = new TH1F("numpadsHist",TString::Format("Number of Hit Pads Vs. Event Number, Run %d", tpcRun), 200, 0, 199);
	  TH1F *charges = new TH1F("chargesHist",TString::Format("Total Event Charge Vs. Event Number, Run %d", tpcRun), 200, 0, 199);
      for (int i = 0; i < 200; ++i) { 
         adc->SetBinContent(i, maxADCs[i]);
		 numpads->SetBinContent(i, totalPadsHit[i]);
		 charges->SetBinContent(i, totalCharges[i]);
      }
	  
      TCanvas *c1 = new TCanvas();
      adc->Draw();
	  
	  TCanvas *c2 = new TCanvas();	  
	  charges->Draw();
	  
	  TCanvas *c3 = new TCanvas();	  
	  numpads->Draw();
      
}