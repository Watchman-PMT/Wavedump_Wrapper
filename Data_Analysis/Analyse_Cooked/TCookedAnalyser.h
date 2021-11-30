#ifndef TCookedAnalyser_h
#define TCookedAnalyser_h

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TRandom3.h>

#include <vector>
#include <limits.h>

#include <numeric>


using namespace std;

class TCookedAnalyser {
 public :

  TFile * inFile = nullptr;

  TFile * outFile = nullptr; 

  // meta data tree for 
  // storing constants
  TTree * metaTree;
  
  short  SampFreq;
  short  NSamples;
  short  NADCBins;
  short  Range_V;
  float  nsPerSamp;
  float  mVPerBin;
  float  Length_ns;
  short  FirstMaskBin;
  float  AmpGain;
  char   FileID[128]; 
  
  int    Run;
  int    PMT;
  int    Loc;
  char   Test;
  int    HVStep;

  TBranch * b_SampFreq     = 0;
  TBranch * b_NSamples     = 0;
  TBranch * b_NADCBins     = 0;
  TBranch * b_Range_V      = 0;
  TBranch * b_nsPerSamp    = 0;
  TBranch * b_mVPerBin     = 0;
  TBranch * b_Length_ns    = 0;
  TBranch * b_AmpGain      = 0;
  TBranch * b_FirstMaskBin = 0;
  TBranch * b_FileID       = 0;
  TBranch * b_Run          = 0;
  TBranch * b_PMT          = 0;
  TBranch * b_Loc          = 0;
  TBranch * b_Test         = 0;
  TBranch * b_HVStep       = 0;

  //--------------------
  // cooked data
  TTree * cookedTree;
  
  // Input 
  vector <short> * ADC = 0;   
  float peak_mV;
  short peak_samp;
  float min_mV;
  float mean_mV;
  float start_s;
  float base_mV;

  TBranch * b_ADC = 0;
  TBranch * b_peak_mV  = 0;  
  TBranch * b_peak_samp = 0;  
  TBranch * b_min_mV  = 0;  
  TBranch * b_mean_mV = 0;  
  TBranch * b_start_s = 0;  
  TBranch * b_base_mV = 0;  
  
  TCookedAnalyser(string path);
  ~TCookedAnalyser();
  int  GetEntry(int entry);
  void Init();

  void  InitMeta();
  void  InitCooked();
   
  void InitCanvas(float w = 1000.,
		  float h = 800.);
  void DeleteCanvas();
  
  string GetFileID();
  
  string GetCookedTreeID();
  string GetMetaTreeID();
  
  void  PrintMetaData();
  
  // limit entries for faster testing
  void  SetTestMode(int);
  
  //---
  // Monitor Noise
  // Noise
  
  float  thresh_mV;
  float  th_low_mV;
  float  noise_thresh_mV;
  float  noise_th_low_mV;
  
  TH1F * hMean_Cooked = nullptr;
  TH1F * hPPV_Cooked  = nullptr;
  
  TH1F * hMin_Cooked = nullptr;
  TH1F * hPeak_Cooked = nullptr;
  
  TH2F * hMin_Peak_Cooked = nullptr;
  
  void  SetTest(char Test);
  char  GetTest();

  void  SetRun(int Run);
  int   GetRun();
  
  void  Noise();
  void  InitNoise();
  void  SaveNoise(string outFolder = "./Plots/Noise/");
  
  float ADC_To_Wave(short ADC);
  float Wave_To_Amp_Scaled_Wave(float wave);

  //----
  void  Make_hQ_Fixed();
  
  bool  HasLowNoise(float min_mV,float peak_mV,
		    float thresh_mV = 10.);

  //----
  // LED data  

  void  Fit_Peak_Time_Dist();
  void  Set_LED_Delay(float LED_delay);
  float Get_LED_Delay();
  
  //----
  // Dark Counts

  TH1F * hD_Peak = nullptr;
  TH2F * hD_Min_Peak = nullptr;

  double base_average(int iEntry);
  double base(int iEntry);
  //TH1D * lin_reg1 = nullptr;
  //TH1D * lin_reg2 = nullptr;
  double average;
  int peak_rise(float thesh_mV = 10., int nbins = 10);
  void  Dark(float thresh_mV = 10.);
  void  InitDark();
  void  SaveDark(string outFolder = "./Plots/Dark/");
  
  void DarkPlot(char option = 'D');

  //---
  // Monitor Waveforms

  TH1F * hWave = nullptr;   
  TH1F * hFFT = nullptr;
  
  TRandom3 * rand3 = nullptr;

  void  Waveform(char option = 'f');
  
  void  InitWaveform();
  void  SaveWaveform(string outPath = "./Plots/Waveforms/");
  
  void  InitFFT();
  void  SaveFFT(string outPath = "./Plots/Waveforms/",
		int    option  = 0);
  
  void  SaveWaveFFT(string outPath = "./Plots/Waveforms/");

  void  SavePulseFit(string outPath = "./Plots/PulseFit/");

  //--- 
  // Rise and Fall Times
   
  TF1 * Fit_Pulse(int entry = -1);
  
  bool IsGoodPulseFit(TF1* f1);
  
/*   float Get_Pulse_Rise_Time(); */
/*   float Get_Pulse_Fall_Time(); */
  
 private:
  
   // only accommodating int size here
   Long64_t nentries64_t; // dummy
   int      nentries;

   TCanvas * canvas = nullptr;
   
   void  Set_THF_Params(float *,float *,float *, int *);
   
   void  SetStyle();
   float fLED_Delay;
   bool IsTimingDistFitted;
   
};

#endif

#ifdef TCookedAnalyser_cxx
TCookedAnalyser::TCookedAnalyser(string path) 
{

  inFile = new TFile(path.c_str(),"READ");
  
  if ( !inFile || !inFile->IsOpen()) {
    fprintf(stderr,"\n Error, Check File: %s \n",path.c_str());
    return;
  }

  Init();
}

TCookedAnalyser::~TCookedAnalyser()
{ 
  delete cookedTree;
  delete metaTree;
  delete inFile;
}

int TCookedAnalyser::GetEntry(int entry)
{
// Read contents of entry.
   if (!cookedTree) return 0;
   return cookedTree->GetEntry(entry);
}

void TCookedAnalyser::SetTestMode(int user_nentries = 1000000){

  nentries = user_nentries;  
  printf("\n Warning: \n ");
  printf("  nentries set to %d for testing \n",nentries);
  
}

void TCookedAnalyser::Init()
{

  InitMeta();
  InitCooked();
  
  rand3 = new TRandom3(0);

  SetStyle();  
  InitCanvas();
  
}

void TCookedAnalyser::InitMeta(){
  
  printf("\n ------------------------------ \n");
  printf("\n Initialising Meta Data \n");
  printf("\n   %s \n",GetMetaTreeID().c_str());

  inFile->GetObject(GetMetaTreeID().c_str(),metaTree);
  
  if (!metaTree){
    fprintf( stderr, "\n Error: no meta tree  \n ");
    fprintf( stderr, "\n Was this file created with the latest cook_raw ? \n ");
    return;
  }

  metaTree->SetBranchAddress("SampFreq",&SampFreq,&b_SampFreq);
  metaTree->SetBranchAddress("NSamples",&NSamples,&b_NSamples);
  metaTree->SetBranchAddress("NADCBins",&NADCBins,&b_NADCBins);
  metaTree->SetBranchAddress("Range_V",&Range_V,&b_Range_V);
  metaTree->SetBranchAddress("nsPerSamp",&nsPerSamp,&b_nsPerSamp);
  metaTree->SetBranchAddress("mVPerBin",&mVPerBin,&b_mVPerBin);
  metaTree->SetBranchAddress("Length_ns",&Length_ns,&b_Length_ns);
  metaTree->SetBranchAddress("AmpGain",&AmpGain,&b_AmpGain);
  metaTree->SetBranchAddress("FirstMaskBin",&FirstMaskBin,&b_FirstMaskBin);
  metaTree->SetBranchAddress("FileID",&FileID,&b_FileID);
  
  metaTree->SetBranchAddress("Run",&Run,&b_Run);
  metaTree->SetBranchAddress("PMT",&PMT,&b_PMT);
  metaTree->SetBranchAddress("Loc",&Loc,&b_Loc);
  metaTree->SetBranchAddress("Test",&Test,&b_Test);
  metaTree->SetBranchAddress("HVStep",&HVStep,&b_HVStep);

  metaTree->GetEntry(0);
  
  printf("\n ------------------------------ \n");
  printf("\n FileID = %s ",FileID);
  printf("\n Run    = %d ",Run);
  printf("\n PMT    = %d ",PMT);
  printf("\n Loc    = %d ",Loc);
  printf("\n Test   = %c ",Test);
  printf("\n HVStep = %d \n",HVStep);

  fLED_Delay = 100.;
  IsTimingDistFitted = false;
  
  printf("\n ------------------------------ \n");

}

void TCookedAnalyser::InitCooked(){
  
  inFile->GetObject(GetCookedTreeID().c_str(),cookedTree);
  
  if (cookedTree == 0){
    fprintf( stderr, "\n Warning: No cooked data tree");
  }
  
  printf("\n ------------------------------ \n");
  printf("\n Initialising Cooked Data \n");
  printf("\n   %s \n",GetCookedTreeID().c_str());
  
  if (!cookedTree){
    fprintf( stderr, "\n Error: no cooked tree  \n ");
    return;
  }

  cookedTree->SetMakeClass(1);
  
  cookedTree->SetBranchAddress("ADC",&ADC, &b_ADC);
  cookedTree->SetBranchAddress("peak_mV",&peak_mV, &b_peak_mV);
  cookedTree->SetBranchAddress("peak_samp",&peak_samp, &b_peak_samp);
  cookedTree->SetBranchAddress("min_mV",&min_mV, &b_min_mV);
  cookedTree->SetBranchAddress("mean_mV",&mean_mV, &b_mean_mV);
  cookedTree->SetBranchAddress("start_s",&start_s, &b_start_s);
  cookedTree->SetBranchAddress("base_mV",&base_mV, &b_base_mV);
  
  nentries64_t = cookedTree->GetEntriesFast();
  
  if( nentries64_t > INT_MAX ){
      fprintf(stderr,
	      "\n Error, nentries = (%lld) > INT_MAX unsupported \n ",
	      nentries64_t);
      return;
  }
  else
    nentries = (int)nentries64_t;
  
  
  printf("\n ------------------------------ \n");
  
  return;
}

void TCookedAnalyser::SetTest(char userTest ){
  Test = userTest;
}

char TCookedAnalyser::GetTest(){
  return Test;
}

void TCookedAnalyser::SetRun(int userRun ){
  Run = userRun;
}

int TCookedAnalyser::GetRun(){
  return Run;
}

#endif // #ifdef TCookedAnalyser_cxx
