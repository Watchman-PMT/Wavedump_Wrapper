#define TCookedAnalyser_cxx
#include "TCookedAnalyser.h"

#include "TSystem.h"
#include <TH2.h>
#include <TF1.h>
#include <math.h>
#include <limits.h>
#include <fstream>

#include "wmStyle.C"

string TCookedAnalyser::GetCookedTreeID(){
  string CookedTreeID = FileID;
  return "Cooked_" + CookedTreeID;
}

string TCookedAnalyser::GetMetaTreeID(){
  return "Meta_Data";
}

string TCookedAnalyser::GetFileID(){
  return FileID;
}

//----------
//   

void TCookedAnalyser::Make_hQ_Fixed(){

  float gate_start = Get_LED_Delay() - 15;

  
  float gate_width = 50.; 
  
  string fileName = "hQ_Fixed_";
  string histName = "hQ_Fixed_";
  
  fileName += GetFileID();
  histName += GetFileID();
  
  fileName += ".root";
  
  printf("\n  %s \n",fileName.c_str());
  
  outFile = new TFile(fileName.c_str(),"RECREATE",fileName.c_str()); 
  
  int   nBins = 100.;
  //float minQ  = -100.0;
  float minQ  = -1000.0;
  float maxQ  = 1000.0;

  TH1F * hQ_Fixed = new TH1F(histName.c_str(),"hQ_Fixed;Charge (mV ns);Counts",
			     nBins,minQ,maxQ);

  float wave_mV = 0.0;
  float time_ns = 0.0;

  int   nSigSamps = 0;
  int   nBasSamps = 0;
  
  float sig_volts = 0;
  float bas_volts = 0;

  float volts  = 0.0;
  float charge = 0.0;

  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    cookedTree->GetEntry(iEntry);
    
    // zero at start of event
    sig_volts = 0.0;
    bas_volts = 0.0;
    nSigSamps = 0;
    nBasSamps = 0;
    
    for( short iSamp = 0 ; iSamp < NSamples; iSamp++){
      
      wave_mV = ADC_To_Wave(ADC->at(iSamp));
      time_ns = iSamp*nsPerSamp;

      ///........
      // in development
      //
      // Sum ADC values in pulse region 
      // and subtrace baseline here. 
      // NB pulse pol has already been made positive.
      if     ( time_ns >=  gate_start   && 
	       time_ns <  (gate_start + gate_width) ){
	sig_volts += wave_mV; 
	nSigSamps++;
      }// baseline subtraction
      else if( time_ns <   gate_start   &&
      	       time_ns >= (gate_start - gate_width ) ){ 
      	bas_volts += wave_mV;
      	nBasSamps++;
      }
      
      volts = sig_volts - (bas_volts*(float)nSigSamps/nBasSamps);
    }
    
    // Convert ADC to units of charge then
    // fill histogram 
    charge = volts*nsPerSamp;

    //if(HasLowNoise(min_mV,peak_mV))
    hQ_Fixed->Fill(charge);
  
  }
  
  hQ_Fixed->Draw();
  gPad->SetLogy();
  
  histName = "./Plots/Charge/" + histName;
  histName += ".pdf";
  
  gPad->SaveAs(histName.c_str());
  
  gPad->SetLogy(false);
  
  outFile->cd();
  outFile->Write();
  outFile->Delete();

}

bool TCookedAnalyser::HasLowNoise(float min_mV, float peak_mV,
				  float thresh_mV){
  
  if     (min_mV < -2.5 && peak_mV < thresh_mV)
    return false;
  else if(peak_mV < -2*min_mV && peak_mV > thresh_mV )
    return false;
  else 
    return true;
}


//
//----------

void TCookedAnalyser::Fit_Peak_Time_Dist(){
  
  printf("\n ------------------------------ \n");
  printf("\n Getting LED delay   \n");

  InitCanvas();
    
  float binWidth = nsPerSamp;
  float hist_min_time = 0.0;
  float hist_max_time = 200.0;  
  float fit_min_time  = 50.0;
  float fit_max_time  = 150.0;
  int   nBins = 0;

  // fix the binning
  Set_THF_Params(&hist_min_time,&hist_max_time,&binWidth,&nBins);
  
  TH1F * hPeakTime = new TH1F("hPeakTime",
			      "hPeakTime;peak time (ns);Counts",
			      nBins,hist_min_time,hist_max_time);
  float LED_delay   = 1.0;
  float delay_width = 1.0;
  float peak_time = 0.;
  float thresh_mV = 10.;
  
  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    cookedTree->GetEntry(iEntry);
    
    peak_time = peak_samp * nsPerSamp;
    
    if(peak_mV > thresh_mV)
      hPeakTime->Fill(peak_time);

  }

  hPeakTime->Fit("gaus","QR","",fit_min_time,fit_max_time);
  
  TF1 * fPeak = hPeakTime->GetFunction("gaus");
  LED_delay   = fPeak->GetParameter(1);
  delay_width = fPeak->GetParameter(2);
  fPeak->Delete();

  fit_min_time = LED_delay - 2*delay_width; 
  fit_max_time = LED_delay + 2*delay_width; 

  hPeakTime->Fit("gaus","QR","",fit_min_time,fit_max_time);
  
  hPeakTime->GetFunction("gaus")->SetLineColor(kBlue);

  hPeakTime->Draw();  

  string outPath = "./Plots/Timing/";
  string sys_command = "mkdir -p ";
  sys_command += outPath;
  gSystem->Exec(sys_command.c_str());

  canvas->SaveAs("./Plots/Timing/hPeakTime.pdf");
 
  DeleteCanvas();

  LED_delay   = hPeakTime->GetFunction("gaus")->GetParameter(1);
  delay_width = hPeakTime->GetFunction("gaus")->GetParameter(2);
  
  printf("\n delay = %.1f (%.1f) \n", LED_delay, delay_width);

  IsTimingDistFitted = true;
  Set_LED_Delay(LED_delay);
}

void TCookedAnalyser::Set_LED_Delay(float LED_delay){
  fLED_Delay = LED_delay;
}

float TCookedAnalyser::Get_LED_Delay(){
  
  if(IsTimingDistFitted)
    return fLED_Delay;
  else {
    Fit_Peak_Time_Dist();
    return Get_LED_Delay();
  }
}

//------------------------------
void TCookedAnalyser::Noise(){
  
  InitNoise();
  
  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    cookedTree->GetEntry(iEntry);
    
    //printf(" \n mean_mV = %f \n", mean_mV);

    hMean_Cooked->Fill(mean_mV);
    hPPV_Cooked->Fill(peak_mV-min_mV);
    hPeak_Cooked->Fill(peak_mV);
    hMin_Cooked->Fill(min_mV);
    hMin_Peak_Cooked->Fill(min_mV,peak_mV);
    
   }// end: for (int iEntry = 0;
  
  // find peak of mean voltage in mV
  int     max_bin_mean = hMean_Cooked->GetMaximumBin();
  TAxis * x_axis       = hMean_Cooked->GetXaxis();
  float   peak_mean_mV = x_axis->GetBinCenter(max_bin_mean);
  
  thresh_mV       = 10.0; // ideally 1/4 of 1 p.e.
  th_low_mV       = 5.0;  // 
  
  noise_thresh_mV = peak_mean_mV - thresh_mV;
  noise_th_low_mV = peak_mean_mV - th_low_mV;
  
  // standard threshold rel mean peak
  int thresh_bin   = hMin_Cooked->FindBin(noise_thresh_mV);
  int noise_counts = hMin_Cooked->Integral(0,thresh_bin);
  
  float noise_rate = (float)noise_counts/nentries;
  noise_rate = noise_rate/Length_ns * 1.0e9;

  // low threshold rel mean peak
  thresh_bin   = hMin_Cooked->FindBin(noise_th_low_mV);
  noise_counts = hMin_Cooked->Integral(0,thresh_bin);
  
  float noise_rate_low = (float)noise_counts/nentries;
  noise_rate_low = noise_rate_low/Length_ns * 1.0e9;
  
  printf("\n Mean voltage %.2f mV \n",peak_mean_mV);
  printf("\n Rate @ %.2f mV  %.2f Hz \n",noise_th_low_mV,noise_rate_low);
  printf("\n Rate @ %.2f mV  %.2f Hz \n",noise_thresh_mV,noise_rate);

  SaveNoise();

}

void TCookedAnalyser::InitNoise(){
  
  printf("\n ------------------------------ \n");
  printf("\n Analysing Noise   \n");
  
  float range = (float)roundf(Range_V)*1000.;
  
  // for amp not 10x, scale bin width 
  float binWidth =  Wave_To_Amp_Scaled_Wave(mVPerBin);
  float minX     = -range/2.;
  float maxX     =  range/2.;
  int   nBins    = 0;
  
  // fix binning and set number of bins
  Set_THF_Params(&minX,&maxX,&binWidth,&nBins);
  
  hMean_Cooked = new TH1F("hMean_Cooked",
			  ";mean voltage (mV);Counts",
			  nBins,minX,maxX);

  hPeak_Cooked =  new TH1F("hPeak_Cooked",
			   ";peak voltage (mV);Counts",
			   nBins,minX,maxX);
  
  hMin_Cooked =  new TH1F("hMin_Cooked",
			  ";min voltage (mV);Counts",
			  nBins,minX,maxX);
  
  //   printf("\n nBins    = %d \n",nBins);
//   printf("\n minX   = %f \n",minX);
//   printf("\n maxX   = %f \n",maxX);
//   printf("\n binWidth = %f \n",binWidth);
//   printf("\n mVPerBin = %f \n",mVPerBin);

  hMin_Peak_Cooked =  new TH2F("hMin_Peak_Cooked",
			       "peak vs min ; min voltage (mV);peak voltage (mV)",
			       nBins,minX,maxX,
			       nBins,minX,maxX);

  // prepare for range starting at zero
  minX = 0.0;
  maxX = range/2.;
  nBins  = 0;
  
  Set_THF_Params(&minX,&maxX,&binWidth,&nBins);
  
  hPPV_Cooked =  new TH1F("hPPV_Cooked",
			  ";peak to peak voltage (mV);Counts",
			  nBins,minX,maxX);

}


void TCookedAnalyser::SaveNoise(string outPath){

  string sys_command = "mkdir -p ";
  sys_command += outPath;
  gSystem->Exec(sys_command.c_str());
  
  printf("\n Saving Noise Monitoring Plots \n\n");

  InitCanvas();
  
  gPad->SetLogy();
  
  hMean_Cooked->SetAxisRange(-30., 120.,"X");
  hMean_Cooked->SetMinimum(0.1);
  hMean_Cooked->Draw();

  string outName = outPath + "hMean_Cooked.pdf";
  canvas->SaveAs(outName.c_str());  
  
  hPPV_Cooked->SetAxisRange(-5.0, 145.,"X");
  hPPV_Cooked->SetMinimum(0.1);
  hPPV_Cooked->Draw();
  
  outName = outPath + "hPPV_Cooked.pdf";
  canvas->SaveAs(outName.c_str());
  
  hPeak_Cooked->SetAxisRange(-20.,80.,"X");
  hPeak_Cooked->SetMinimum(0.1);
  hPeak_Cooked->Draw();
  outName = outPath + "hPeak_Cooked.pdf";
  canvas->SaveAs(outName.c_str());
  
  
  hMin_Cooked->SetAxisRange(-30.,20.,"X");
  hMin_Cooked->SetMinimum(0.1);
  hMin_Cooked->Draw();

  TLine * l_thresh = new TLine(noise_thresh_mV,1,noise_thresh_mV,1000);
  l_thresh->SetLineStyle(2);
  l_thresh->SetLineColor(kRed);
  l_thresh->SetLineWidth(2);
  l_thresh->Draw();
  
  TLine * l_th_low = new TLine(noise_th_low_mV,1,noise_th_low_mV,1000);
  l_th_low->SetLineStyle(2);
  l_th_low->SetLineColor(kBlue);
  l_th_low->SetLineWidth(2);
  l_th_low->Draw();
  
  outName = outPath + "hMin_Cooked.pdf";
  canvas->SaveAs(outName.c_str());
  
  gPad->SetLogy(false);
  gPad->SetLogz();
  
  hMin_Peak_Cooked->SetAxisRange(-25., 15.,"X");
  hMin_Peak_Cooked->SetAxisRange(-15., 50.,"Y");
  
  hMin_Peak_Cooked->Draw("colz");

  outName = outPath + "hMin_Peak_Cooked.pdf";
  canvas->SaveAs(outName.c_str());

  gPad->SetLogz(false);
  
  DeleteCanvas();
  
}

double TCookedAnalyser::base_average(int iEntry){

  std::vector<double> amplitude;

  for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
    amplitude.push_back(ADC_To_Wave(ADC->at(iSamp)));
    
  double sum = std::accumulate(amplitude.begin(), amplitude.end(), 0.0);
  double mean = sum/NSamples;
  
  return mean;

}

double TCookedAnalyser::base(int iEntry){

  std::vector<double> amplitude;

  for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
    amplitude.push_back(ADC_To_Wave(ADC->at(iSamp)));
    
  double sum = std::accumulate(amplitude.begin(), amplitude.end(), 0.0);
  double mean = sum/NSamples;
  
  //need baseline subtraction here
  
  std::vector<double> base_sub;
  
  for( int isamp = 0; isamp < NSamples; isamp++){
    base_sub.push_back(amplitude[isamp]-mean);
  }
  
  double var = 0.;
  
  for (int n = 0; n < NSamples; n++){
    var += pow(base_sub[n]-mean,2);
  }

  var /= NSamples;
  double sdev = sqrt(var);
  
  std::vector<double> av_base;
  int nSamp = 0;
  for (int N = 0; N < NSamples; N++){
    if (sqrt(pow(base_sub[N]-mean,2)) <= 3*sdev){
      av_base.push_back(base_sub[N]);
      nSamp++;
    }
  }
  
  double baseline = std::accumulate(av_base.begin(), av_base.end(), 0.0);
  baseline /= nSamp;

  return baseline;

}

int TCookedAnalyser::peak_rise(float thresh_mV, int nbins){

  double thresh = base_mV+0.5*peak_mV;//base_mV + thresh_mV;
  
  std::vector<double> amplitude;

  for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
    amplitude.push_back(ADC_To_Wave(ADC->at(iSamp)));
    
  int bins = 0;
    
  for( int iSamp_peak = peak_samp; iSamp_peak > peak_samp - nbins; iSamp_peak--){
    if(amplitude[iSamp_peak] > thresh)
      bins++;
    else
      break;
  }
    
  if(bins == 0)
    return 0;  
  else if(bins == nbins)
    return 0;
  else
    return 1;
    
  //first analysis uses 6 bins, base_mV + thresh_mV, 10 mV thresh, no bins == 0 condition
  
}

void TCookedAnalyser::Dark(float thresh_mV){
  
  InitDark();
  
  float darkRate = 0;
  float darkRateErr = 0;
  float darkRate_noise = 0;
  float darkRateErr_noise = 0;
  
  TFile* results = new TFile("dark_results.root","RECREATE");  
  TTree* Dark = new TTree("Dark","Dark");
  Dark->Branch("darkRate",&darkRate,"darkRate/F");
  Dark->Branch("darkRateErr",&darkRateErr,"darkRateErr/F");
  Dark->Branch("darkRate_noise",&darkRate_noise,"darkRate/F");
  Dark->Branch("darkRateErr_noise",&darkRateErr_noise,"darkRateErr/F");
  
  int nDark = 0;
  int nDark_noise = 0;
  
  int rejected = 0;
  
  int rise_rej = 0;
  int av_neg_rej = 0;
  int av_pos_rej = 0;
  int peak_low = 0;
  int peak_high = 0;
  
  std::ofstream rejected_waveforms;
  rejected_waveforms.open("rejected_waveforms.csv");
  rejected_waveforms << "Rejected waveform at entry\n";

  std::ofstream dark_csv;
  dark_csv.open ("dark_hits.csv");
  dark_csv << "Count at entry\n";
  
  for (int iEntry = 0; iEntry < nentries; iEntry++) {
    cookedTree->GetEntry(iEntry);
     
    if(peak_mV > thresh_mV)
      nDark_noise++;

    // Noise Rejection 
    if( min_mV < -2.5 && peak_mV < thresh_mV)
      continue;
    
    if( peak_mV < -2*min_mV && peak_mV > thresh_mV )
      continue;
      
    if( peak_mV < 2*min_mV && peak_mV > thresh_mV )
      continue;
      
    hD_Peak->Fill(peak_mV);
    hD_Min_Peak->Fill(min_mV,peak_mV);
    
    if( peak_mV < thresh_mV){
      peak_low++;
      continue;}
    
    average = base_average(iEntry);
    
    if( average < -10){
      rejected_waveforms << iEntry << "\n";
      rejected++;
      av_neg_rej++;
      continue;}
      
    if( average > 10){
      rejected_waveforms << iEntry << "\n";
      rejected++;
      av_pos_rej++;
      continue;}
      
    if( peak_mV > 100){
      rejected_waveforms << iEntry << "\n";
      rejected++;
      peak_high++;
      continue;}
    
    int rise = peak_rise();
    
    if(!rise){
      rejected_waveforms << iEntry << "\n";
      //rejected++;
      rise_rej++;
      continue;}
    
    dark_csv << iEntry << "\n";
    
    nDark++;
    
  }

  rejected_waveforms.close();
  dark_csv.close();
  
  std::ofstream rej_count;
  rej_count.open("rejected_types.csv");
  rej_count << "peak_low,av_neg_rej,av_pos_rej,peak_high,rise_rej\n";
  rej_count << peak_low << "," << av_neg_rej << "," << av_pos_rej << "," << peak_high << "," << rise_rej;
  rej_count.close();

  float darkErr = sqrt(nDark);

  darkRate = (float)nDark/(nentries-rejected);
  darkRate = darkRate/Length_ns * 1.0e9;
  darkRateErr = darkErr/nDark * darkRate;
  
  printf("\n \n nentries = %d \n",nentries);
  printf("\n %i rejected 'dark counts'\n",rejected);
  printf("\n dark counts (noise rejected) = %d +/- %.0f \n",nDark,darkErr);
  printf("\n dark rate   (noise rejected) = %.0f +/- %.0f Hz \n",darkRate,darkRateErr);
  
  std::ofstream dark_results;
  dark_results.open ("dark_results.txt");
  dark_results << "dark counts (noise rejected) = " << nDark << " +/- " << darkErr << "\n";
  dark_results << "dark noise (noise rejected) = " << darkRate << " +/- " << darkRateErr << " Hz\n";
  dark_results.close();
  
  float darkErr_noise = sqrt(nDark_noise);
  
  darkRate_noise = (float)nDark_noise/nentries;
  darkRate_noise = darkRate_noise/Length_ns * 1.0e9;
  darkRateErr_noise = darkErr_noise/nDark_noise * darkRate_noise;
  
  printf("\n dark counts (with noise) = %d +/- %.0f \n",nDark_noise,darkErr_noise);
  printf("\n dark rate   (with noise) = %.0f +/- %.0f Hz\n\n",darkRate_noise,darkRateErr_noise);
  
  std::ofstream dark_results_noise;
  dark_results_noise.open ("dark_results.txt", std::ios_base::app);
  dark_results_noise << "dark counts (noise) = " << nDark_noise << " +/- " << darkErr_noise << "\n";
  dark_results_noise << "dark noise (noise) = " << darkRate_noise << " +/- " << darkRateErr_noise << " Hz\n";
  dark_results_noise.close();
  
  SaveDark();
  
  Dark->Fill();
  Dark->Write();
  results->Close();
  
}

void TCookedAnalyser::InitDark(){
  
  printf("\n ------------------------------ \n");
  printf("\n Dark Counts Analysis           \n");
    
  float range = (float)roundf(Range_V)*1000.;

  float max      =  range/2;
  float min      = -range/2;
  float binWidth = Wave_To_Amp_Scaled_Wave(mVPerBin);
  int   nBins    = 0;

  //  fix binning and set number of bins
  Set_THF_Params(&min,&max,&binWidth,&nBins);
  
//   hBase = new TH1F("hBase",
// 		   "hBase;baseline voltage (mV);Counts",
// 		   nBins,min,max);
  
  hD_Peak = new TH1F("hD_Peak",
		     "hD_Peak;peak voltage (mV);Counts",
		     nBins,min,max);
  
//   hBase_Peak = new TH2F("hBase_Peak",
// 			"hBase_Peak;baseline voltage (mV);peak voltage (mV)",
// 			nBins,min,max,
// 			nBins,min,max);
  
  hD_Min_Peak = new TH2F("hD_Min_Peak",
			 "hD_Min_Peak;min voltage (mV);peak voltage (mV)",
			 nBins,min,max,
			 nBins,min,max);

}


void TCookedAnalyser::SaveDark(string outPath){

  string sys_command = "mkdir -p ";
  sys_command += outPath;
  gSystem->Exec(sys_command.c_str());

  InitCanvas();

  TLegend *leg = new TLegend(0.21,0.2,0.31,0.9);
    
  leg->SetTextSize(0.025);
  leg->SetHeader("Baseline start","C");
  
  leg->SetMargin(0.4); 

  gPad->SetLogy();
  
  hD_Peak->SetAxisRange(-5., 75.,"X");
  hD_Peak->SetMinimum(0.1);
  hD_Peak->Draw();
  
  TLine * lVert = new TLine(10,0,10,20);
  lVert->SetLineColor(kBlue);
  lVert->SetLineWidth(2);
  lVert->SetLineStyle(2);
  lVert->Draw();

  string outName = outPath + "hD_Peak.pdf";
  canvas->SaveAs(outName.c_str());

  gPad->SetLogy(false);

//   hBase_Peak->SetAxisRange(-25.,25.,"X");
//   hBase_Peak->SetAxisRange(-5., 45.,"Y");
  
//   hBase_Peak->Draw("col");
  
  
//   outName = outPath + "hBase_Peak.pdf";
//   canvas->SaveAs(outName.c_str());

  gPad->SetLogz();
  hD_Min_Peak->SetAxisRange(-25.,25.,"X");
  hD_Min_Peak->SetAxisRange(-5., 45.,"Y");
  
  gPad->SetGrid(1, 1);
  hD_Min_Peak->Draw("col");
  
  gPad->SetLogz();
  
  outName = outPath + "hD_Min_Peak.pdf";
  canvas->SaveAs(outName.c_str());

  gPad->SetGrid(0,0);
  
  gPad->SetLogz(false);

  DeleteCanvas();
  
}

float TCookedAnalyser::ADC_To_Wave(short ADC){

  float wave = ADC * mVPerBin;

  wave -= Range_V*1000./2.;
  
  wave = Wave_To_Amp_Scaled_Wave(wave);
  
  return wave;
}

float TCookedAnalyser::Wave_To_Amp_Scaled_Wave(float wave){
  return wave/AmpGain*10.;
}


void TCookedAnalyser::DarkPlot(char option){

  //fix FFT, add option to cycle through dark and rejected plots rather than plotting all

  char answer = 'D';
  
  printf("\n What to plot? \n");
  printf("\n D - All dark waveforms \n");
  printf("\n R - All rejected dark waveforms \n");
  printf("\n X - eXit \n");
  // note deliberate use of whitespace before %c
  scanf(" %c", &answer);
  
  char plot = 'b';
  
  printf("\n How to save plots? \n");
  printf("\n p - png \n");
  printf("\n r - root \n");
  printf("\n b - 'p' and 'r' \n");
  scanf(" %c", &plot);

  if(answer == 'D'){

    if (plot!='p'){
      string fileName = "hWaveDark_";
      fileName += GetFileID();
      fileName += ".root";
      printf("\n  %s \n",fileName.c_str());
      outFile = new TFile(fileName.c_str(),"RECREATE",fileName.c_str());
    }
  
    //InitFFT();
    InitWaveform();
    
    std::ifstream file ("dark_hits.csv");
    std::string line;
    std::getline(file,line); //remove first line
    std::vector<int> indices;
        
    while(std::getline(file,line))
    {
      std::string line_value;   
      std::stringstream ss(line);
      while(std::getline(ss,line_value,','))
      {
        indices.push_back(stoi(line_value));
        
      }
    }
           
    int nWaveforms = indices.size();
    
    for (int i = 0; i < nWaveforms ; i++){
      cookedTree->GetEntry(indices[i]);
      
      for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
        hWave->SetBinContent(iSamp+1,(ADC_To_Wave(ADC->at(iSamp))));
      
      if(plot == 'p'){
        
        string outPath = "./Plots/Waveforms/Dark/Kept/";
        string sys_command = "mkdir -p ";
        sys_command += outPath;
        gSystem->Exec(sys_command.c_str());
      
        string outPathfile = outPath;
        outPathfile += "hWave_";
        outPathfile += to_string(indices[i]);
        outPathfile += ".png";
        SaveWaveform(outPathfile);
      }
      else if(plot == 'r'){
        outFile->Write();
      }
      else{
        
        string outPath = "./Plots/Waveforms/Dark/Kept/";
        string sys_command = "mkdir -p ";
        sys_command += outPath;
        gSystem->Exec(sys_command.c_str());
        
        string outPathfile = outPath;
        outPathfile += "hWave_";
        outPathfile += to_string(indices[i]);
        outPathfile += ".png";
        SaveWaveform(outPathfile);
        
        outFile->Write();
      }
    }
    
  }
  
  if(answer == 'R'){
  
    if(plot!='p'){
      string fileName = "hWaveDarkRejected_";
      fileName += GetFileID();
      fileName += ".root";
      printf("\n  %s \n",fileName.c_str());
      outFile = new TFile(fileName.c_str(),"RECREATE",fileName.c_str());
    }
    
    //InitFFT();
    InitWaveform();
    
    std::ifstream file ("rejected_waveforms.csv");
    std::string line;
    std::getline(file,line); //remove first line
    std::vector<int> indices;
        
    while(std::getline(file,line))
    {
      std::string line_value;   
      std::stringstream ss(line);
      while(std::getline(ss,line_value,','))
      {
        indices.push_back(stoi(line_value));
        
      }
    }
    
    int nWaveforms = indices.size();
    
    for (int i = 0; i < nWaveforms ; i++){
      cookedTree->GetEntry(indices[i]);
      
      for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
        hWave->SetBinContent(iSamp+1,(ADC_To_Wave(ADC->at(iSamp))));
          
      if(plot == 'p'){
      
        string outPath = "./Plots/Waveforms/Dark/Rejected/";
        string sys_command = "mkdir -p ";
        sys_command += outPath;
        gSystem->Exec(sys_command.c_str());
      
        string outPathfile = outPath;
        outPathfile += "hWave_";
        outPathfile += to_string(indices[i]);
        outPathfile += ".png";
        SaveWaveform(outPathfile);
      }
      else if(plot == 'r'){
        outFile->Write();
      }
      else{
      
        string outPath = "./Plots/Waveforms/Dark/Rejected/";
        string sys_command = "mkdir -p ";
        sys_command += outPath;
        gSystem->Exec(sys_command.c_str());
      
        string outPathfile = outPath;
        outPathfile += "hWave_";
        outPathfile += to_string(indices[i]);
        outPathfile += ".png";
        SaveWaveform(outPathfile);
      
        outFile->Write();
      }
    }
    
  }

}

void TCookedAnalyser::Waveform(char option){

  string fileName = "hWave_";
  fileName += GetFileID();
  fileName += ".root";
  printf("\n  %s \n",fileName.c_str());
  outFile = new TFile(fileName.c_str(),"RECREATE",fileName.c_str()); 
  
  switch(option){
  case('w'):
    InitFFT();// inc Waveform
    break;
  case('f'):
    InitFFT();// inc Waveform
    break;
  case('b'):
    InitFFT();// inc Waveform
    break;
  default:
    //InitWaveform();
    InitFFT();
    break;
  }
  
  int  entry  = 0;
  char answer = 'R';
  
  // Plotting Loop
  while ( answer!='X' ){
    
    if(option!='d'){
      printf("\n Which entry to plot? \n");
      printf("\n F - First entry \n");
      printf("\n N - Next \n");
      printf("\n P - Previous \n");
      printf("\n R - Random selection \n");
      printf("\n A - Accumulate All \n");
      printf("\n I - Index of waveform \n");
      printf("\n X - eXit \n");
      
      // note deliberate use of whitespace before %c
      scanf(" %c", &answer);
    }
    
    switch(answer){
    case('F'):
      entry = 0;
      break;
    case('R'):
      entry = (int)round(rand3->Uniform(nentries)); ;
      break;
    case('N'):
      entry++;
      break;
    case('P'):
      entry--;
      break;
    case('A'):
      entry = 0;
      break;
    case('I'):
      std::cout << "Entry to plot: \n" << endl;
      std::cin >> entry;
      break;
    default:
      entry = -1;
    }
    
    if(entry > -1)
      
      if( answer=='A' )
	printf("\n averaging all entries \n");
      else
	printf("\n plotting entry %d \n",entry);
    
    else{
      printf("\n exiting waveform plotting \n");
      return;
    }
    
    cookedTree->GetEntry(entry);
    
    // single or initial instance
    for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
      hWave->SetBinContent(iSamp+1,(ADC_To_Wave(ADC->at(iSamp))));
    
    hWave->FFT(hFFT ,"MAG");
    
    // accumulate
    int nWaveforms = 1;
    
    if(answer == 'A' )
      nWaveforms  = nentries;
    
    TH1F * hWave_temp = nullptr; 
    TH1F * hFFT_temp  = nullptr;

    for (int i = (entry+1); i < nWaveforms ; i++){
      cookedTree->GetEntry(i);

      if( i%100000 == 0 )
	printf("\n entry %d/%d \n",i,nWaveforms );
      
      hWave_temp = (TH1F*)hWave->Clone();
      hFFT_temp  = (TH1F*)hFFT->Clone();
      
      for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
	hWave_temp->SetBinContent(iSamp+1,(ADC_To_Wave(ADC->at(iSamp))));
      
      hWave_temp->FFT(hFFT_temp ,"MAG");
      
      hWave->Add(hWave_temp);
      hFFT->Add(hFFT_temp);
      
      hWave_temp->Delete();
      hFFT_temp->Delete();
    }

    hFFT->SetBinContent(1,0.);

    double_t norm = 1;
    hFFT->Scale(norm/hFFT->Integral(),"width");

    hFFT->SetMaximum(0.1);
    hFFT->SetMinimum(0.0001);

    //outFile->cd();
    outFile->Write();
    //outFile->Delete();

    string outPath = "./Plots/Waveforms/";

    string sys_command = "mkdir -p ";
    sys_command += outPath;
    gSystem->Exec(sys_command.c_str());

    switch(option){
    case('w'):
      outPath += "hWave.pdf";
      SaveWaveform(outPath);
      break;
    case('f'):
      outPath += "hFFT.pdf";
      SaveFFT(outPath);
      //SaveFFT(outPath,1); // SetLogx() option
      break;
    case('b'):
      outPath += "hWaveFFT.pdf";
      SaveWaveFFT(outPath);
      break;
    case('d'):
      // make default plots
      outPath = "./Plots/Waveforms/";
      outPath += "hWaveFFT.pdf";
      SaveWaveFFT(outPath);
      
      outPath = "./Plots/Waveforms/";
      outPath += "hWave.pdf";
      SaveWaveform(outPath);
      
      outPath = "./Plots/Waveforms/";
      outPath += "hFFT.pdf";
      SaveFFT(outPath);
      //SaveFFT(outPath,1); // SetLogx() option    
     
      answer='X';
      break;
    default:
      break;
    }
    
  }
  
}

void TCookedAnalyser::InitWaveform(){
  
  printf("\n ------------------------------ \n");
  printf("\n Init Wave \n\n");
  
  hWave = new TH1F("hWave","Waveform;Time (ns); Amplitude (mV)",
		   NSamples, 0.,Length_ns);
  
}

void TCookedAnalyser::InitFFT(){

  printf("\n ------------------------------ \n");
  printf("\n Init FFT \n\n");
  
  hWave = new TH1F("hWave","Waveform;Time (ns); Amplitude (mV)",
		   NSamples, 0.,Length_ns);

  hFFT = new TH1F("hFFT","FFT; Frequency (MHz); Magnitude",
		  NSamples/2, 0, SampFreq/2 );
  
}


void TCookedAnalyser::SaveWaveform(string outPath ){


  printf("\n Saving Waveform Plot \n\n");
  
  InitCanvas();
  
  hWave->Draw();
  
  canvas->SaveAs(outPath.c_str());

  //hWave->Delete();
  //  hFFT->Delete();
  DeleteCanvas();
  
}

void TCookedAnalyser::SaveFFT(string outPath, int option){
  
  printf("\n Saving FFT Plot \n\n");
  
  InitCanvas();
  
  
  hFFT->SetAxisRange(5,450.,"X");
  gPad->SetLogx(option);
  //gPad->SetLogy(option);
  hFFT->Draw("hist");
  
  canvas->SaveAs(outPath.c_str());
  
  //  hWave->Delete();
  //hFFT->Delete();
  DeleteCanvas();
  
}

void TCookedAnalyser::SaveWaveFFT(string outPath){
  
  printf("\n Saving Waveform and FFT Plots \n\n");
  
  InitCanvas(1600.);
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hWave->Draw();
  
  canvas->cd(2);
  hFFT->SetBinContent(1,0.);
  hFFT->Draw("hist");
  
  canvas->SaveAs(outPath.c_str());

  //hWave->Delete();
  //hFFT->Delete();
  DeleteCanvas();
  
}


TF1 * TCookedAnalyser::Fit_Pulse(int entry){
  
  float LED_delay = Get_LED_Delay(); 
  float min_time_range = LED_delay - 65;  // 50
  float max_time_range = LED_delay + 185; // 300
  
  if(max_time_range > Length_ns)
    max_time_range = Length_ns;
  
  if(min_time_range < 0)
    min_time_range = 0.0;
  
  // "crystalball" params defined below (fNumber = 500)
  // https://root.cern.ch/root/html608/TFormula_8cxx_source.html
  // https://en.wikipedia.org/wiki/Crystal_Ball_function
  // 0 - const (height of peak)
  // 1 - mean  
  // 2 - sigma (width)
  // 3 - alpha (small appears to give bigger skew)
  // 4 - N     (normalisation)
  
  // add baseline patameter
  string fName = "[0]+crystalball(1)";
  
  TF1 * fWave = new TF1("fWave",fName.c_str(),
			min_time_range,max_time_range);
  
  int attempts = 0;
  do{
    
    attempts++;
    
    //if(entry==-1)
    entry = (int)round(rand3->Uniform(nentries));

    printf("\n nentries =  %d \n",nentries);
    printf("\n Fitting entry %d \n",entry);
    
    cookedTree->GetEntry(entry);
    
    InitWaveform();

    for( short iSamp = 0 ; iSamp < NSamples; iSamp++)
      hWave->SetBinContent(iSamp+1,(ADC_To_Wave(ADC->at(iSamp))));
    
    // basline,height,mean,sigma,alpha,norm
    fWave->SetParameters(base_mV,10,LED_delay,2.5,-1,1);
    
    fWave->SetParLimits(0, base_mV-5,base_mV+5); // baseline
    fWave->SetParLimits(1, 0, 1000); // height
    fWave->SetParLimits(2, LED_delay-7.5,LED_delay+7.5);//mean
    fWave->SetParLimits(3,0.5,5); // sigma
    fWave->SetParLimits(4,-2,-0.5); // alpha
    fWave->FixParameter(5,1); // N
    
    hWave->SetAxisRange(min_time_range,max_time_range,"X");
    
    hWave->Fit("fWave", "R");
  }
  while( !IsGoodPulseFit(fWave) );
  
  hWave->Draw();

  printf("\n There were %d attempts to fit. \n",attempts);
  
  SavePulseFit();
  
  return fWave;
}

bool TCookedAnalyser::IsGoodPulseFit(TF1 * f1){
  
  //  float base   = f1->GetParameter(0);
  float height = f1->GetParameter(1);
//   float mean   = f1->GetParameter(2);
//   float sigma  = f1->GetParameter(3);
//   float norm   = f1->GetParameter(4); 
  
  float xi_sq  = f1->GetChisquare()/f1->GetNDF(); 
  
  if ( height < 10 ){
    printf("\n pulse amplitude is %f mV \n",height);
    return false;
  }
  else if( xi_sq > 5 || xi_sq < 0.01 ){
    printf("\n #chi^{2} is %f  \n",xi_sq);
    return false;
  }
  else{
    printf("\n #chi^{2}/NDF is %f  \n",xi_sq);
    return true;
  }
}

void TCookedAnalyser::SavePulseFit(string outPath ){
  
  string sys_command = "mkdir -p ";
  sys_command += outPath;
  gSystem->Exec(sys_command.c_str());

  printf("\n Saving Pulse Fit Plot \n\n");
  
  InitCanvas();
  
  hWave->Draw();
  TF1 * f = hWave->GetFunction("fWave");
  f->Draw("same");

  outPath += "Pulse.pdf";

  canvas->SaveAs(outPath.c_str());

  DeleteCanvas();
  
}



//------------------------------
//Fix histogram binning
void TCookedAnalyser::Set_THF_Params(float * minX, 
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

void TCookedAnalyser::InitCanvas(float w,float h){
  
  canvas = new TCanvas();
  canvas->SetWindowSize(w,h);
  
}

void TCookedAnalyser::DeleteCanvas(){
  delete canvas;
}

void TCookedAnalyser::PrintMetaData(){ 

  printf("\n ------------------------------ \n");
  printf("\n Printing Meta Data \n");

  if (!metaTree) return;
  metaTree->Show(0);
  
}

void TCookedAnalyser::SetStyle(){
  
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
