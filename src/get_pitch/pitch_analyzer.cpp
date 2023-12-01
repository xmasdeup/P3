/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"
#include <vector>
#include <complex>

using namespace std;

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
  
  // Función para calcular la FFT recursivamente
  void PitchAnalyzer::calcularFFT(std::vector<std::complex<double>> &x, bool inverse) const{
    const int N = x.size();
    double exponent = 0;
    if (N <= 1) return;

    std::vector<std::complex<double>> even(N/2), odd(N/2);
    for (int i = 0; i < N/2; ++i) {
        even[i] = x[2*i];
        odd[i] = x[2*i + 1];
    }

    calcularFFT(even,inverse);
    calcularFFT(odd,inverse);
    if(inverse == false)  exponent =  -2 * M_PI / N;
    else if(inverse == true) exponent =  2 * M_PI / N;

    for (int k = 0; k < N/2; ++k) {
        std::complex<double> t = std::polar(1.0,k*exponent) * odd[k];
        x[k] = even[k] + t;
        x[k + N/2] = even[k] - t;
    }
  }

  // Función para calcular el cepstrum
  void PitchAnalyzer::calcularCepstrum(const std::vector<float> &signal, std::vector<float> &cepstrum) const {
    int N = cepstrum.size();
    // Transformada de Fourier de la señal
    std::vector<std::complex<double>> fftSignal(N);

    for (int i = 0; i < N; ++i) {
        fftSignal[i] = { signal[i], 0.0 }; // Se llena con valores complejos (la parte imaginaria es cero)
    }
    
    calcularFFT(fftSignal,false);

    // Calculando el logaritmo del espectro
    std::vector<double> logMagnitude(N);
    for (int i = 0; i < N; ++i) {
        logMagnitude[i] = std::log(std::abs(fftSignal[i]));
    }

    // Transformada inversa de Fourier del logaritmo del espectro
    std::vector<std::complex<double>> ifftLogMagnitude(N);

    std::copy(logMagnitude.begin(), logMagnitude.end(), ifftLogMagnitude.begin());
    calcularFFT(ifftLogMagnitude, true); // FFT inversa
    float_t maxIndex = 0;
    int index = 0;
    // Extrayendo el cepstrum (parte real de la IFFT del logaritmo del espectro)
    for (int i = 0; i < N; ++i) {
        cepstrum[i] = ifftLogMagnitude[i].real();
        if(maxIndex<cepstrum[i]) 
        {
          maxIndex = cepstrum[i];
          index = i;
        }
    }
    // std::cout<<maxIndex<<" \t";
    // std::cout<<index<<" \n";
  }

  std::float_t PitchAnalyzer::compute_zcr(const std::vector<float> &x) const
  {
    std::float_t count = 0;

    for(uint64_t i = 1; i < x.size(); i++)
    {
        if(x[i]<x[i-1])
        {
          count++;
        }
    }

    return count;
  }

  void PitchAnalyzer::amdf(const vector<float> &x, vector<float> &distance) const{
      
      float_t minIndex = 99999;
      float_t mean = 0;
      int index = 0;
      std::cout<<x.size()<<"\n";
      for (unsigned int lag = npitch_min; lag < distance.size() ;lag++)
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
      std::cout<<mean<<"\t";
      std::cout<<minIndex<<" \t";
      std::cout<<index<<" \n";
      std::cout<<"End Of Frame\n";

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

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    if((pot< -40) || (rmaxnorm/r1norm < 0.35))
    {
      return true;
    } 
    else return false;
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);
    vector<float> distance(npitch_max);
    vector<float> c(npitch_max);
    
    //Compute correlation
    autocorrelation(x, r);
    calcularCepstrum(x, c);
    amdf(x,distance);
    vector<float>::const_iterator iR = r.begin(), iRMax = iR;
    
    for(iR = iRMax = r.begin() + npitch_min; iR< r.begin() + npitch_max; iR++)
    {
      if(*iR>*iRMax)
      {
        iRMax = iR;
      }
    }

    unsigned int lag = iRMax - r.begin();

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
  ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.


    float pot = 10 * log10(r[0]);
    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif 
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0]))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }

  ButterWorthFilter::ButterWorthFilter(int32_t order, double_t cutoffFrequency, double_t samplingFrequency):
      order(order), cutoffFrequency(cutoffFrequency), samplingFrequency(samplingFrequency) {

        double_t wc = 2.0 * M_PI * cutoffFrequency;

        poles.resize(order);

        for(int32_t k = 0; k<order; ++k)
        {
          double_t theta = (2 * k + 1) * M_PI / (2.0*order);
          double_t realPart = -sin(theta) * sinh(1.0/(2.0*order)*asinh(1.0));
          double_t imagPart = cos(theta) * cosh(1.0/(2.0*order)*asinh(1.0)); 
          poles[k] =  std::complex<double>(realPart,imagPart) * wc;
        } 

        for(int32_t k = 0; k<order;++k)
        {
          poles[k]= (2.0 * samplingFrequency + poles[k]) / (2.0 * samplingFrequency - poles[k]);
        }

        gain = 1.0;
        for(int32_t k= 0; k< order;++k)
        {
          gain *= -poles[k].real();
        }

  }
  void ButterWorthFilter::applyFilter(std::vector<float> &spectrum)
  {
    int N = spectrum.size() /2;

    for(int32_t i = 0; i< N ;++i) {

      double_t freq = static_cast<double>(i) * samplingFrequency/N;

      std::complex<double> H = 1.0;

      for(int32_t k = 0; k<order;++k)
      {
        double poleReal = poles[k].real();
        double poleImag = poles[k].imag();
        std::complex<double> term = 1.0 - exp(-2.0 * M_PI * freq / samplingFrequency * std::complex<double>(poleReal,poleImag));
        H /= term;
      }
      H /= gain;

      spectrum[i] *= H.real();
      spectrum[N+i] *= H.imag();
    }
  }

  void ButterWorthFilter::center_clipping(std::vector<float> &signal) const
  {
    float_t maximum = 0;
    for(int32_t i = 0; i<signal.size();++i)
    {
      if(abs(signal[i]) > maximum) maximum = abs(signal[i]);
    }

    float_t threshold = maximum * 0.3;

    for(int32_t k = 0; k<signal.size();++k)
    {
      if(signal[k] >= (threshold)) signal[k] -= threshold;
      else if(signal[k] <= (-threshold)) signal[k] += threshold;
      else signal[k] = 0;
    }
    

  }

}
