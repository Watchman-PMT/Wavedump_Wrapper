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

#include "FileNameParser.h"

using namespace std;

class TCookedAnalyser {
 public :

  FileNameParser * fnp = nullptr;
  TFile * inFile = nullptr;

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
  char   FileID_char[128]; 
  string FileID; 
  
  TBranch * b_SampFreq  = 0;
  TBranch * b_NSamples  = 0;
  TBranch * b_NADCBins  = 0;
  TBranch * b_Range_V   = 0;
  TBranch * b_nsPerSamp = 0;
  TBranch * b_mVPerBin  = 0;
  TBranch * b_Length_ns = 0;
  //TBranch * b_FileID    = 0;

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
  
  char  GetTest();
  
  void  Noise();
  void  InitNoise();
  void  SaveNoise(string outFolder = "./Plots/Noise/");
  
  float ADC_To_Wave(short ADC);

  //----
  // LED data
  
  float Get_LED_delay();

  //----
  // Dark Counts

  TH1F * hD_Peak = nullptr;
  TH2F * hD_Min_Peak = nullptr;

  void  Dark(float thresh_mV = 10.);
  void  InitDark();
  void  SaveDark(string outFolder = "./Plots/Dark/");

  //---
  // Monitor Waveforms

  TH1F * hWave = nullptr;   
  TH1F * hFFT = nullptr;
  
  TRandom3 * rand3 = nullptr;

  void  Waveform(char option = 'f');
  
  void  InitWaveform();
  void  SaveWaveform(string outFolder = "./Plots/Waveforms/");
  
  void  InitFFT();
  void  SaveFFT(string outFolder = "./Plots/Waveforms/");
  
  void  SaveWaveFFT(string outFolder = "./Plots/Waveforms/");

 private:
  
   // only accommodating int size here
   Long64_t nentries64_t; // dummy
   int      nentries;

   TCanvas * canvas = nullptr;
   
   void  Set_THF_Params(float *,float *,float *, int *);
   
   void  SetStyle();

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
  
  // set FileID using inFile path
  fnp = new FileNameParser(path);  
  FileID = fnp->GetFileID();
  
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
    return;
  }

  metaTree->SetBranchAddress("SampFreq",&SampFreq,&b_SampFreq);
  metaTree->SetBranchAddress("NSamples",&NSamples,&b_NSamples);
  metaTree->SetBranchAddress("NADCBins",&NADCBins,&b_NADCBins);
  metaTree->SetBranchAddress("Range_V",&Range_V,&b_Range_V);
  metaTree->SetBranchAddress("nsPerSamp",&nsPerSamp,&b_nsPerSamp);
  metaTree->SetBranchAddress("mVPerBin",&mVPerBin,&b_mVPerBin);
  metaTree->SetBranchAddress("Length_ns",&Length_ns,&b_Length_ns);
  //metaTree->SetBranchAddress("FileID",FileID_char,&b_FileID);
  
  //sprintf(FileID,"%s",FileID_char);

  metaTree->GetEntry(0);
  

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

char TCookedAnalyser::GetTest(){
  
  if(fnp)
    return fnp->GetTest();
  else
    return 'E';
}


#endif // #ifdef TCookedAnalyser_cxx
