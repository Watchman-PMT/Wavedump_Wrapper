#define TCooker_cxx
#include "TCooker.h"
#include <TH2.h>
#include <math.h>
#include <limits.h>

#include "wmStyle.C"

void TCooker::Cook(){
  
  // initialise trees
  InitCooking();

  // create variables in standard units
  // and find waveform peak
  DoCooking();

  SaveMetaData();
  SaveCookedData();
  outFile->Close();
    
  printf("\n Cooking is complete            \n");
  printf("\n ------------------------------   ");
  printf("\n ------------------------------ \n");
  
}

void TCooker::InitCooking(){

  printf("\n ------------------------------ \n");
  printf("\n Initialising Cook \n");

  InitCookedDataFile();
  InitMetaDataTree();
  InitCookedDataTree();

}


void TCooker::InitCookedDataFile(string option){

  string fileName = GetDir();
  fileName += GetFileID();
  fileName += ".root";
    
  if     (!strcmp(option.c_str(),"RECREATE")) 
    printf("\n Preparing: ");
  else if(!strcmp(option.c_str(),"OPEN")) 
    printf("\n Opening:   ");
  else   {
    printf("\n Error: option unknown \n");
    return;
  }
    
  printf("\n  %s \n",fileName.c_str());

  outFile = new TFile(fileName.c_str(),
		      option.c_str(),
		      fileName.c_str());

}

void TCooker::SaveMetaData(){
  
  printf("\n ------------------------------ \n");
  printf("\n Writing meta data              \n");   
  
  sprintf(FileID,"%s",f_fileID.c_str());

  metaTree->Fill();

  metaTree->Write();
  metaTree->Delete();
  
}

void TCooker::SaveCookedData(){
  
  printf("\n ------------------------------ \n");
  printf("\n Writing cooked data            \n");
  printf("\n Closing:                         ");
  printf("\n   %s       \n",outFile->GetName());
  printf("\n ------------------------------ \n");
    
  cookedTree->Write();
  cookedTree->Delete();


}

void TCooker::InitCookedDataTree(){
  
  // ----------
  // Temporary tree for calibrating
  string treeName = "Cooked_";
  treeName += GetFileID();

  cookedTree = new TTree(treeName.c_str(),treeName.c_str());

  cookedTree->Branch("ADC",&ADC_buff);
  cookedTree->Branch("peak_mV",&peak_mV,"peak_mV/F");
  cookedTree->Branch("peak_samp",&peak_samp,"peak_samp/S");
  cookedTree->Branch("min_mV",&min_mV,"min_mV/F");
  cookedTree->Branch("mean_mV",&mean_mV,"mean_mV/F");
  cookedTree->Branch("start_s",&start_s,"start_s/F");
  cookedTree->Branch("base_mV",&base_mV,"base_mV/F"); 
  
}

void TCooker::InitMetaDataTree(){
  
  string treeName = "Meta_Data";

  //string treeName = "Meta_";
  //treeName += GetFileID();

  metaTree = new TTree(treeName.c_str(),treeName.c_str());
  
  metaTree->Branch("SampFreq",&fSampFreq,"SampFreq/S");
  metaTree->Branch("NSamples",&fNSamples,"NSamples/S");
  metaTree->Branch("NADCBins",&fNADCBins,"NADCBins/S");
  metaTree->Branch("Range_V",&fRange_V,"Range_V/S");
  metaTree->Branch("nsPerSamp",&f_nsPerSamp,"nsPerSamp/F");
  metaTree->Branch("mVPerBin",&f_mVPerBin,"mVPerBin/F");
  metaTree->Branch("Length_ns",&fLength_ns,"Length_ns/F");
  metaTree->Branch("AmpGain",&fAmpGain,"AmpGain/F");
  metaTree->Branch("FirstMaskBin",&fFirstMaskBin,"FirstMaskBin/S");  
  metaTree->Branch("FileID",FileID,"FileID/C");

  //
  metaTree->Branch("Run",&fRun,"Run/I");  
  metaTree->Branch("PMT",&fPMT,"PMT/I");  
  metaTree->Branch("Loc",&fLoc,"Loc/I");  
  metaTree->Branch("Test",&fTest,"Test/B");  
  metaTree->Branch("HVStep",&fHVStep,"HVStep/I");  
  //
}

void TCooker::DoCooking(){
  
  printf("\n ------------------------------ \n");
  printf("\n Cooking                       \n");
  
  int    nBaseSamps;
  double time = 0, prevTime = 0; 
  int    trigCycles = 0;
    
  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    rawTree->GetEntry(iEntry);
  
    wave_mV.clear(), ADC_buff.clear();
    nBaseSamps = 0,     peak_samp  =  0   ;
    min_mV     = 1000., peak_mV   = -1000.;
    base_mV    = 0.,    mean_mV   =  0.   ;

    // event start time
    time = GetElapsedTime(&trigCycles,prevTime);
    prevTime = time; // now set for next entry
    start_s = (float)time; 
    
    // first loop - find baseline, set wave_mV
    for (short iSamp = 0; iSamp < fNSamples; ++iSamp){
      // voltage scaled to pre-amp gain
      // plus pulse flip if necessary
      wave_mV.push_back(ADC_To_Wave(ADC->at(iSamp)));      
      if( IsSampleInBaseline(iSamp) ){
	base_mV += wave_mV.at(iSamp);
	nBaseSamps++;
      }
    }
    base_mV /= (float)nBaseSamps;
    
    // second loop - apply baseline subtraction, set variables
    for (short iSamp = 0; iSamp < fNSamples; ++iSamp){
      
      // subtract baseline or mask 
      if( fFirstMaskBin > 0 && 
	  iSamp > fFirstMaskBin){
	wave_mV.at(iSamp) = 0.0;
      }
      else{
	wave_mV.at(iSamp) -= base_mV;
      }
      // min
      if( wave_mV.at(iSamp) < min_mV)
	min_mV = wave_mV.at(iSamp);
      // peak_samp
      if( wave_mV.at(iSamp) >= peak_mV ){
	peak_mV = wave_mV.at(iSamp);
	peak_samp  = iSamp;
      }
      mean_mV += wave_mV.at(iSamp);
      
      // ADC
      // ADC mask
      if( fFirstMaskBin > 0 && 
	  iSamp >= fFirstMaskBin  ){ 
	ADC_buff.push_back(Wave_To_ADC(base_mV));	
      } // ADC flip
      else{ 
	// if pulse polarity is negative then flip 
	ADC_buff.push_back(Invert_Negative_ADC_Pulses(ADC->at(iSamp)));
      }
    }
    mean_mV = mean_mV/(float)fNSamples;    
    
    cookedTree->Fill();
  }
  
}


// void TCooker::SetFileID(){
//   f_fileID = "fileID";
// }

void TCooker::SetFileID(string userFileID){
  
  f_fileID = userFileID;
  
  printf("\n TCooker object FileID set to: ");
  printf("\n  %s \n ",f_fileID.c_str());
  
}


void TCooker::SetDir(string userFileDir){
  
  f_fileDir = userFileDir;
  
  printf("\n TCooker object FileDir set to: ");
  printf("\n  %s \n ",f_fileDir.c_str());
  
}

//
void TCooker::SetRun(int userRun){
  
  fRun = userRun;
  
  printf("\n TCooker object Run set to: ");
  printf("\n  %d \n ",fRun);
  
}

void TCooker::SetPMT(int userPMT){
  fPMT = userPMT;
}

void TCooker::SetLoc(int userLoc){
  fLoc = userLoc;
}

void TCooker::SetTest(char userTest){
  fTest = userTest;
}

void TCooker::SetHVStep(int userHVStep){
  fHVStep = userHVStep;
}

string TCooker::GetFileID(){
  return f_fileID;
}

string TCooker::GetDir(){
  return f_fileDir;
}

string TCooker::GetCookedTreeID(){
  return "Cooked_" + GetFileID();  
}

short TCooker::Invert_Negative_ADC_Pulses(short ADC){
  
  if(fPulsePol=='N'){
    ADC -= GetNADCBins()/2;
    ADC = -ADC; 
    ADC += GetNADCBins()/2;
  }
  
  return ADC;
}

float TCooker::ADC_To_Wave(short ADC){

  float wave = ADC * Get_mVPerBin();

  wave -= GetRange_mV()/2.;

  if(fPulsePol=='N')
    wave = -wave;

  wave = Wave_To_Amp_Scaled_Wave(wave);
  
  return wave;
}

short TCooker::Wave_To_ADC(float wave){

  wave = Amp_Scaled_Wave_To_Unscaled_Wave(wave);

  if(fPulsePol=='N')
    wave = -wave;

  wave += GetRange_mV()/2.;
  
  wave = wave/Get_mVPerBin();
  
  short ADC = (short)roundf(wave);

  ADC = Invert_Negative_ADC_Pulses(ADC);
  
  return ADC; 
}


float TCooker::Wave_To_Amp_Scaled_Wave(float wave){
  return wave/fAmpGain*10.;
}

float TCooker::Amp_Scaled_Wave_To_Unscaled_Wave(float wave){
  return wave*fAmpGain/10.;
}

bool TCooker::IsSampleInBaseline(short iSample){
  // adjust this function to use most populous voltage bin(s)?
  // use ADC_To_Wave to convert ADC to mV
  float sampTime = (float)iSample * SampleToTime(); //convert sample to time
  float width    = 50.; // sets baseline to first 50 ns
  
  if( sampTime <= width)
    return true;
  else
    return false;
}

/*float TCooker::Baseline(short ADC){
// need ADC in
// find most populous bin
// convert this to voltage
// set this as baseline
// create bool function like IsSampleInBaseline (?)

}*/

//------------------------------

void TCooker::InitBaseline(){
  
  printf("\n ------------------------------ \n");
  printf("\n Analysing Baseline           \n\n");
  
  float mVPerBin = Get_mVPerBin();
  float low_mV   = -GetRange_mV()/2.;
  float high_mV   =  GetRange_mV()/2.;
  int   nBins    = 0;
  
  // fix binning and set number of bins
  Set_THF_Params(&low_mV,&high_mV,&mVPerBin,&nBins);
  
  hBase = new TH1F("hBase",
		   "hBase;baseline voltage (mV);Counts",
		   nBins,low_mV,high_mV);
  
  hPeak = new TH1F("hPeak",
		   "hPeak;peak voltage (mV);Counts",
		   nBins,low_mV,high_mV);
  
  hBase_Peak = new TH2F("hBase_Peak",
			"hBase_Peak;baseline voltage (mV);peak voltage (mV)",
			nBins,low_mV,high_mV,
			nBins,low_mV,high_mV);
  
  hMin_Peak = new TH2F("hMin_Peak",
		       "hMin_Peak;min voltage (mV);peak voltage (mV)",
		       nBins,low_mV,high_mV,
		       nBins,low_mV,high_mV);

  float minEvent = 0;
  float maxEvent = float(nEvents_Base-1);
  float eventsPerBin = 0.; // dummy
  int   nEvents;

  minEvent = 0.;
  maxEvent = 0. + (float)nEvents_Base - 1.;
  eventsPerBin = 0.;
  nEvents = nEvents_Base;
    
  Set_THF_Params(&minEvent,&maxEvent,&eventsPerBin,&nEvents);

  hEvent_Base = new TH2F("hEvent_Base",
			 ";Event;baseline voltage (mV)",
			 nEvents,minEvent,maxEvent,
			 nBins,low_mV,high_mV);
 
}

void TCooker::SaveBaseline(string outFolder){

  printf("\n Saving Baseline Study Plots \n\n");

  InitCanvas();

  TLegend *leg = new TLegend(0.21,0.2,0.31,0.9);
    
  //  char title[128] = "";
  leg->SetTextSize(0.025);
  leg->SetHeader("Baseline start","C");
  
  //   Int_t colors[] = {kRed-7,kRed,kRed+2,
  //     kRed-5,kOrange,kOrange+2,
  //     kOrange+4,kYellow+1,kGreen+1,
  //     kGreen+3,kGreen-5,kCyan+1,
  //     kCyan+3,kBlue-7,kBlue,
  //     kBlue+3,kViolet,kMagenta+1};
  
  leg->SetMargin(0.4); 


  hBase->SetAxisRange(-25., 25.,"X");
  hBase->SetMinimum(0.1);

  hBase->Draw();
  
  //  sprintf(title,"%.0f ns",iBase*10.);
  //  leg->AddEntry(hBase[iBase],title,"L");  
  //  leg->Draw();
  
  gPad->SetLogy();
  string outName = outFolder + "hBase.pdf";
  canvas->SaveAs(outName.c_str());
  
  hPeak->SetAxisRange(-5., 75.,"X");
  hPeak->SetMinimum(0.1);
  hPeak->Draw();

  outName = outFolder + "hPeak.pdf";
  canvas->SaveAs(outName.c_str());
  
  gPad->SetLogy(false);

  hBase_Peak->SetAxisRange(-25.,25.,"X");
  hBase_Peak->SetAxisRange(-5., 65.,"Y");
  
  hBase_Peak->Draw("col");
  
  gPad->SetLogz();
 
  outName = outFolder + "hBase_Peak.pdf";
  canvas->SaveAs(outName.c_str());
  
  hMin_Peak->SetAxisRange(-25.,5.,"X");
  hMin_Peak->SetAxisRange(-5., 65.,"Y");
  
  gPad->SetGrid(1, 1);
  hMin_Peak->Draw("col");
  
  TLine * lVert = new TLine(-2.5,0,-2.5,10.);
  TLine * lDiag = new TLine(-5,10,-20,40.);
  
  lVert->SetLineStyle(2);
  lVert->SetLineColor(kBlue);
  lVert->SetLineWidth(2);

  lDiag->SetLineStyle(2);
  lDiag->SetLineColor(kBlue);
  lDiag->SetLineWidth(2);

  gPad->SetLogz();
    
  lVert->Draw();
  lDiag->Draw();
  
  outName = outFolder + "hMin_Peak.pdf";
  canvas->SaveAs(outName.c_str());

  gPad->SetGrid(0,0);
  
  gPad->SetLogz(false);

  float w = 10000., h = 100.;
  canvas->SetWindowSize(w,h);
  
  gPad->SetGrid(1, 1);
  gPad->SetLogz();
  
  float base_mean = 0.0;
  float minX = -1.5, maxX = 1.5;
  base_mean = hBase->GetMean(); 
  minX = base_mean - 1.5;
  maxX = base_mean + 1.5;
  
  hEvent_Base->SetAxisRange(minX,maxX,"Y");
  hEvent_Base->Draw("col");

  outName = outFolder + "hEvent_Base.pdf";
  canvas->SaveAs(outName.c_str());

  gPad->SetGrid(0, 0);
  w = 1000.;
  h = 800.;
  canvas->SetWindowSize(w,h);
  gPad->SetLogz(false);

  DeleteCanvas();
}

//------------------------------
// Fix histogram binning
void TCooker::Set_THF_Params(float * minX, 
			      float * maxX,
			      float * binWidth,
			      int   * nBins){
  
  if     (*nBins==0)
    *nBins = (int)roundf((*maxX - *minX)/(*binWidth));
  else if(*nBins > 0 && *binWidth < 1.0E-10)
    *binWidth = (*maxX - *minX)/(*nBins);
  else
    fprintf(stderr,"\n Error in Set_THF_Params \n");
  
  *nBins += 1;
  *minX -= 0.5*(*binWidth);
  *maxX += 0.5*(*binWidth);

};

void TCooker::InitCanvas(float w,float h){

  canvas = new TCanvas();
  canvas->SetWindowSize(w,h);

}

void TCooker::DeleteCanvas(){
  delete canvas;
}

double TCooker::GetTrigTimeTag() {

  unsigned int TTT;
   
  // limit to 31 bits 
  if( HEAD[5] > (unsigned int)INT_MAX )
    TTT = HEAD[5] - INT_MAX;
  else
    TTT = HEAD[5];
  
  return (double)TTT*8.E-9;
}

double TCooker::GetTrigTimeTag(int entry) {

  b_HEAD->GetEntry(entry);

  return GetTrigTimeTag();
}

double TCooker::GetElapsedTime(int * cycles,
			       double prevTime) {

  double time = GetTrigTimeTag();
  
  // range is 31 bit with 8 ns LSB = 17.18 seconds 
  // using half range of INT here
  double range = ((double)(INT_MAX)+1.0)*8.E-9;
  
  time += range*(double)(*cycles); 
  
  if(time < prevTime){
    time += range; 
    (*cycles)+= 1;
  }
  
  time -= startTime;
  
  //printf("\n time = %f \n",time);
  
  return time;
}

void TCooker::CountMissedEvents(int dTrigEntry){

  if( dTrigEntry != 1 ){
    nMissedEvents += (dTrigEntry-1);
    printf("\n %d missed events \n",nMissedEvents);
  }
  
}

float TCooker::GetLength_ns(){
  return fLength_ns;
}

short TCooker::GetNSamples(){
  return fNSamples;
}

// private
void TCooker::SetDigitiser(char digitiser){

  if(digitiser=='V' || 
     digitiser=='D')
    fDigitiser = digitiser;
  else{
    fprintf( stderr, "\n Error: unknown digitiser \n ");
    fprintf( stderr, "\n Setting to default ('V')  \n ");
    fDigitiser = 'V';
  }
  
  return;
}


void TCooker::SetAmpGain(float amp_gain){  
  fAmpGain = amp_gain ;
}

void TCooker::SetFirstMaskBin(short first_mask_bin){  
  fFirstMaskBin = first_mask_bin;
}


float TCooker::GetRange_mV(){
  return (float)fRange_V*1000.;
}

float TCooker::Get_mVPerBin(){
  return f_mVPerBin;
}

short TCooker::GetNADCBins(){
  return fNADCBins;
}

void TCooker::SetSampSet(char sampSet){
  
  if  (fDigitiser=='V')
    fSampSet = 'V';
  else
    fSampSet = sampSet;
  
  return;
}

void TCooker::SetPulsePol(char pulsePol){
  
  if(pulsePol == 'N' || 
     pulsePol == 'P')
    fPulsePol = pulsePol;
  else{
    fprintf( stderr, "\n Error: unknown pulse polarity \n ");
    fprintf( stderr, "\n Setting to default ('N')  \n ");
    fPulsePol = 'N';
  }
  
  return;
}

void TCooker::SetConstants(){
  
  printf("\n Setting Constants \n");
  
  fSampFreq = SetSampleFreq();
  fNSamples = SetNSamples();
  fNADCBins = SetNADCBins();
  fRange_V  = SetRange_V();

  // dependent on above
  f_nsPerSamp = Set_nsPerSamp();
  f_mVPerBin  = Set_mVPerBin();
  fLength_ns  = SetLength_ns();

}

void TCooker::PrintConstants(){ 
  
  printf("\n ------------------------------ \n");
  printf("\n  Acquisition Settings          \n");

  if(fDigitiser=='D'){
    printf("\n  desktop digitiser \n");   
    printf("\n  sampling setting     = %c     ",fSampSet);
  }

  printf("\n  sampling frequency   = %d MHz \n",fSampFreq);
  printf("   period per sample   = %.1f ns  \n",f_nsPerSamp);
  printf("\n  samples per waveform = %d     \n",fNSamples);
  printf("   waveform duration   = %.1f ns  \n",fLength_ns);
  printf("\n  number of ADC Bins   = %hd    \n",fNADCBins);
  printf("  ADC range            = %d V     \n",fRange_V);
  printf("   ADC bin width       = %.2f mV  \n",f_mVPerBin);
  printf("\n  pulse polarity       = %c     \n",fPulsePol);
  printf("\n  pre-amp gain         = %.1f   \n",fAmpGain);
  
  if(fFirstMaskBin!=-1)
    printf("\n  first mask bin      = %hd    \n",fFirstMaskBin);
  
  printf(" \n " );
  
}


//------------------------------
// 

 void TCooker::DAQ()
{
  
  nMissedEvents = 0;

  InitDAQ();
  
  double time     = 0;
  double prevTime = 0;
  
  int trigCycles    = 0;
  int trigEntry     = 0;
  int prevTrigEntry = 0;
  int dTrigEntry    = 0;
  
  int nRateEvents = (int)round(nentries/100); // how many to average

  int     deltaEvents = 0;
  double  deltaT      = 0.;
  double  dTime       = 0.; // time between events

  double  eventRate = 0.;
  double  meanRate  = 0.;

  int nbytes = 0, nb = 0;

  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    nb = rawTree->GetEntry(iEntry);   nbytes += nb;    
    
    //-----------------------------
    // Process Header Information
    time = GetElapsedTime(&trigCycles,prevTime);
    
//     printf("\n prevTime = %f \n\n",prevTime);
//     printf("\n time     = %f \n\n",time);
//     printf("\n dTime    = %f \n\n",dTime);
    
    dTime    = time - prevTime; 
    prevTime = time; // now set for next entry
    
    deltaT   += dTime;   // integrated time
    trigEntry = HEAD[4]; // absolute entry number
    
    hTT_EC->Fill(GetTrigTimeTag(),trigEntry);

    deltaEvents++; // integrated event count
    dTrigEntry    = trigEntry - prevTrigEntry;
    prevTrigEntry = trigEntry;
    
    // skip first entry for intergrated
    // and differential variable plots
    if( iEntry==0 ) 
      continue;
    
    hTrigFreq->Fill(1./dTime/1000.);  

    meanRate += 1./dTime/1000.;
    
    CountMissedEvents(dTrigEntry); // any unwritten events?

    // process after event integration period
    if(iEntry%nRateEvents == 0 ){
      
      hNEventsTime->Fill(time/60.,iEntry);
      
      eventRate = deltaEvents/deltaT;
      
      hEventRate->Fill(time/60.,eventRate/1000.);  
    
      // reset integrated variables
      deltaEvents = 0;
      deltaT = 0.;
      
    } // end of: if(iEntry%nRateEvents
    
    // Process Header Information
    //-----------------------------
    
  } // end of: for (int iEntry = 0...

  meanRate = meanRate/nentries;
  
  //printf("\n mean event rate = %.2f kHz \n",meanRate);

  printf("\n Mean trigger frequency is %.2f kHz \n\n",hTrigFreq->GetMean());
  
  SaveDAQ();

  printf("\n ------------------------------ \n");
}

void TCooker::InitDAQ(){
  
  printf("\n ------------------------------ \n");
  printf("\n Getting DAQ Info             \n\n");

  //----
  float minTime    = 0.0;
  float maxTime    = 16.0; // minutes
  float timePerBin = 1./6000.; // 100th of a second per bin
  int   nTimeBins  = 0;
    
  Set_THF_Params(&minTime,&maxTime,&timePerBin,&nTimeBins);
  
  hNEventsTime = new TH1F("hNEventsTime",
			  "hNEventsTime;Time (mins);Event",
			  nTimeBins,minTime,maxTime);
  
  hEventRate = new TH1F("hEventRate",
			"hEventRate;Time (mins);Mean Rate (kHz)",
			nTimeBins,minTime,maxTime);

  //----
  float minFreq    = 0.0; 
  float maxFreq    = 100.0;    // kHz
  float freqPerBin = 1./10000; // 1/10 Hz
  int   nFreqBins  = 0;
  
  Set_THF_Params(&minFreq,&maxFreq,&freqPerBin,&nFreqBins);
  
  hTrigFreq = new TH1F("hTrigFreq",
		       "hTrigFreq; Rate (kHz); Counts",
		       nFreqBins,minFreq,maxFreq);

  //----
  // hTT_EC
  float minClock   = 0.0;
  float maxClock   = ((float)(INT_MAX)+1.0)*8.E-9;
  int   nClockBins = (int)round(UINT_MAX/1000000) + 1; 
  float secsPerClockBin = 0.;
  
  Set_THF_Params(&minClock,&maxClock,&secsPerClockBin,&nClockBins);
  
  rawTree->GetEntry(0);
  float firstEntry = HEAD[4];
  
  rawTree->GetEntry(nentries-1);
  float lastEntry  = HEAD[4];

  float entriesPerBin = 1000.;
  int   nEntryBins    = 0;

  Set_THF_Params(&firstEntry,&lastEntry,&entriesPerBin,&nEntryBins);

  hTT_EC = new TH2F("hTT_EC","hTT_EC;Trigger Time Tag (secs);Entry",
		    nClockBins,minClock,maxClock,
		    nEntryBins,firstEntry,lastEntry);

}

void TCooker::SaveDAQ(string outFolder){
  
  InitCanvas();
  
  int maxBin = hNEventsTime->GetMaximumBin();
  int minBin = 0;

  //maxTime = hNEventsTime->GetXaxis()->GetBinCenter(maxBin);
  maxBin++;
  
  hNEventsTime->GetXaxis()->SetRange(minBin,maxBin);
  hEventRate->GetXaxis()->SetRange(minBin,maxBin);

  hNEventsTime->Draw("HIST P");

  string outName = outFolder + "hNEventsTime.pdf";
  
  canvas->SaveAs(outName.c_str());
  
  maxBin = hTrigFreq->GetMaximumBin();
  float meanFreq_kHz = hTrigFreq->GetMean();
  float minFreq_kHz  = meanFreq_kHz - 2.;
  float maxFreq_kHz  = meanFreq_kHz + 2.;
  
  hTrigFreq->SetAxisRange(minFreq_kHz, maxFreq_kHz,"X");

  hTrigFreq->Draw("hist");

  outName = outFolder + "hTrigFreq.pdf";
  
  canvas->SaveAs(outName.c_str());

  hEventRate->SetMinimum(minFreq_kHz);
  hEventRate->SetMaximum(maxFreq_kHz);
  
  hEventRate->Draw("HIST P");

  outName = outFolder + "hEventRate.pdf";
  canvas->SaveAs(outName.c_str());
  
  hTT_EC->Draw("colz");

  outName = outFolder + "hTT_EC.pdf";
  canvas->SaveAs(outName.c_str());

  DeleteCanvas();
}

short TCooker::SetSampleFreq(){
  
  if(fDigitiser=='D'){
    switch(fSampSet){
    case '0':
	return 5000;
    case '1':
	return 2500;
    case '2':
	return 1000;
    case '3':
	return 750;
    default:
      return 1000;
    }
  }
  else
    return 500; // 'V'
}

float TCooker::Set_nsPerSamp(){
  return 1000./fSampFreq;
}

float TCooker::SampleToTime(){
  return f_nsPerSamp;
}

short TCooker::SetNSamples(){
  
  rawTree->GetEntry(0);   
  
  uint hdrByts = 24;
  uint smpByts = HEAD[0] - hdrByts;

  if(fDigitiser=='V')
    return smpByts/2; // shorts
  else
    return smpByts/4;   // ints
}

float TCooker::SetLength_ns(){
  return f_nsPerSamp*fNSamples;
}

short TCooker::SetNADCBins(){
  
  if( fDigitiser == 'D' )
    return 4096;
  else   
    return 16384;
}

short TCooker::SetRange_V(){

  if(fDigitiser=='V')
    return 2;
  else
    return 1;
  
}

float TCooker::Set_mVPerBin(){
  
  return 1000.*fRange_V/fNADCBins;
  
}



void TCooker::SetStyle(){
  
  printf("\n Setting Style \n");

  TStyle *wmStyle = GetwmStyle();
 
  const int NCont = 255;
  const int NRGBs = 5;
  
  // Color scheme for 2D plotting with a better defined scale 
  double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };          
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
  wmStyle->SetNumberContours(NCont);
 
  gROOT->SetStyle("wmStyle");
  gROOT->ForceStyle();
 
}
