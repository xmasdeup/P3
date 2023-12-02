/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"
#include <vector>
#include <complex>
#include "FFTReal.h"

using namespace std;
using namespace ffft;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l]
      r[l] = 0;
      for (unsigned int m = 0; m < (x.size()-l); ++m){
        r[l] += x[m]*x[m+l];
      }
      r[l] = r[l]/x.size();
    }
      /** \FET Realitzada la funció d'autocorrelació
       * Implementem       
       * - Inicialitzem...
       * - Acumulem...
       * - Dividim...
       * */ 
    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  void PitchAnalyzer::calcularCepstrum(const std::vector<float> &signal, std::vector<float> &cepstrum, int size_fft, FFTReal <float> &fft_first, FFTReal <float> &fft_second) const
  {
    int N = cepstrum.size();

    float in_fourier[size_fft];
    float normal_fourier[size_fft];
    float log_fourier[size_fft/2];
    float out_cepstrum[size_fft/2];

    for(int i = 0; i < size_fft; i++)
    {
      in_fourier[i] = 0;

      if(i < N)
      {
        in_fourier[i] = signal[i];
      }

    }

    fft_first.do_fft(normal_fourier, in_fourier);

    for(int i = 0; i<(size_fft/2); i++)
    {
      log_fourier[i] = log( pow(normal_fourier[i],2) + pow(normal_fourier[i+size_fft/2],2));
    }

    fft_second.do_fft(out_cepstrum,log_fourier);

    for(int i = 0; i<(size_fft/4); i++)
    {
      cepstrum[i] = pow(out_cepstrum[i],2) + pow(out_cepstrum[i+size_fft/4],2);
      cout<<cepstrum[i]<<"\n";
    }

  }



  int PitchAnalyzer::compute_zcr(const std::vector<float> &x) const
  {
    int count = 0;

    for(uint64_t i = 1; i < x.size(); i++)
    {
        if(signbit(x[i]*x[i-1]))
        {
          count++;
        }
    }
    //cout<<count<<"\n";
    return count;
  }

  unsigned int PitchAnalyzer::amdf(const vector<float> &x, vector<float> &distance) const{
      
      float_t minIndex = 1000;
      float_t mean = 0;
      unsigned int index = 0;
      //std::cout<<x.size()<<"\n";
      for (unsigned int lag = npitch_min-10; lag <= distance.size() ;lag++)
      {
          distance[lag] = 0;

        for(unsigned int n = 0; n < (x.size() -lag); ++n)
        {
          distance[lag] += abs(x[n]-x[n+lag]);
        }
        distance[lag] = distance[lag]/(x.size()-lag);

        if(lag<x.size()/2) mean+=distance[lag];
        // std::cout<<lag<<" \t";
        // std::cout<<distance[lag]<<"\n";

        if(lag<x.size()/2) ;
        if(minIndex>distance[lag]) 
        {
          minIndex = distance[lag];
          index = lag;
        }

      }
      if(minIndex < 1.5e-4) index = 0; 

      // std::cout<<mean<<"\t";
      // std::cout<<minIndex<<" \t";
      // std::cout<<index<<" \n";
      // std::cout<<"End Of Frame\n";
    return index;
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
        /// \FET Implement the Hamming window
        window.assign(frameLen, 1);
        for(uint32_t i = 0; i < window.size(); i++)
        {
          window[i] = 1 * (0.53836 - 0.46164 * cos(2 * M_PI * i / (window.size() - 1)));
        }
      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, unsigned int min, float zcr) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    // if((pot< -40) || (rmaxnorm/r1norm < 0.35))
    // {
    //   return true;
    // } 

    if(((min>=npitch_min) && (min<=npitch_max)) && zcr < 60)
    {

      return false;
    }

    else return true;
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x, FFTReal <float> &fft_first, FFTReal <float> &fft_second) const {
    if (x.size() != frameLen)
      return -1.0F;
  int size_fft = 0;

  for(unsigned int size_f = 1; size_f < frameLen; size_f *= 2) size_fft = size_f;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);
    vector<float> distance(npitch_max);
    vector<float> c(size_fft/4);
    
    //Compute correlation
    int zcr = compute_zcr(x);

    autocorrelation(x, r);
    calcularCepstrum(x, c, size_fft, fft_first, fft_second);
    unsigned int min = amdf(x,distance);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;
    
    for(iR = iRMax = r.begin() + npitch_min; iR< r.begin() + npitch_max; iR++)
    {
      if(*iR>*iRMax)
      {
        iRMax = iR;
      }
      // cout<<*iR<<"\n";
    }

    unsigned int lag = iRMax - r.begin();

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
  ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.

    // cout<<npitch_max<<"\n";
    // cout<<npitch_min<<"\n";
    float pot = 10 * log10(1e-8*r[0]);
    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif 
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0], min, zcr))
      return 0;
    else
      return (float) samplingFreq/(float) min;
  }


  void ButterWorthFilter::applyFilter(float spectrum [], const float filter [])
  {
    for(int i = 0; i< samples; i++)
    {
      spectrum[i] *= filter[i];
      spectrum[i+samples] *= filter[i];

    }
  }


  void ButterWorthFilter::center_clipping(std::vector<float> &signal, int rate) const
  {
    float_t interval = rate/20;
    vector<float>::iterator iX;
    float_t maximum = 0;
    float_t threshold = 0;


    for (iX = signal.begin(); iX + interval < signal.end(); iX = iX + interval) 
    {
      maximum = 0;

      for(int32_t i = 0; i < interval; i++)
      {
        if(abs(*(iX+i)) > maximum) maximum = abs(*(iX+i));
      }

      threshold = maximum * 0.25;
      //cout<<maximum<<"\n";
      for(int32_t k = 0; k < interval; k++)
      {
        if( *(iX+k) >= (threshold)) ;//*(iX+k) -= threshold;
        else if(*(iX+k) <= (-threshold)) ; // *(iX+k) += threshold;
        else 
        {
          *(iX+k) = 0;
          //cout<<"ZEROOO\n";
        }
      }

    }
    if(iX < signal.end()){
      int missing = signal.end() - iX;
      for(int32_t k = 0; k < missing; k++)
      {
      if( *(iX+k) >= (threshold)) ;//*(iX+k) -= threshold;
      else if(*(iX+k) <= (-threshold)) ; // *(iX+k) += threshold;
      else 
      {
        *(iX+k) = 0;
        //cout<<"ZEROOO\n";
      }
      }
    }
  }
}
