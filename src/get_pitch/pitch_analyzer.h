/// @file

#ifndef PITCH_ANALYZER_H
#define PITCH_ANALYZER_H

#include <vector>
#include <algorithm>
#include <complex>
#include "FFTReal.h"

#define MAX_FRAME 2048

using namespace ffft;

namespace upc {
  const float MIN_F0 = 50.0F;    ///< Minimum value of pitch in Hertzs
  const float MAX_F0 = 500.0F; ///< Maximum value of pitch in Hertzs

  ///
  /// PitchAnalyzer: class that computes the pitch (in Hz) from a signal frame.
  /// No pre-processing or post-processing has been included
  ///

  class ButterWorthFilter {
  private:

    int samples;
    

  public: 

    ButterWorthFilter(int N){
      samples = N;
    }

  ///
  /// Applies the ButterWorthFilter
  ///
    void applyFilter(float spectrum [], const float filter[]);

  ///
  /// Returns the signal center clipped
  ///
    void center_clipping(std::vector<float> &signal, int rate) const;

  };

  class PitchAnalyzer {
  public:
	/// Wndow type
    enum Window {
		RECT, 						///< Rectangular window
		HAMMING				///< Hamming window

  };

    void set_window(Window type); ///< pre-compute window

  private:
    std::vector<float> window; ///< precomputed window
    unsigned int frameLen, ///< length of frame (in samples). Has to be set in the constructor call
      samplingFreq, ///< sampling rate (in samples per second). Has to be set in the constructor call
      npitch_min, ///< minimum value of pitch period, in samples
      npitch_max, 
      size_fft; ///< maximum value of pitch period, in samples;
    unsigned int previous_min;
    std::vector<float> pitch;

    int fft_frame = MAX_FRAME;

 
	///
	/// Computes correlation from lag=0 to r.size()
	///
    void autocorrelation(const std::vector<float> &x, std::vector<float> &r) const;

  ///
  /// Computes the average maximum distance from lag=1 to MAX_F0 or npitch_max
  ///
    unsigned int amdf(const std::vector<float> &x, std::vector<float> &distance) const;

  ///
  /// Computes the fft of the signal
  ///
    void calcularFFT(std::vector<std::complex<double>> &x, bool inverse) const;


  ///
  /// Computes the cepstrum of the signal
  ///
    void calcularCepstrum(const std::vector<float> &signal, std::vector<float> &cepstrum, int size_fft, FFTReal <float> &fft_first, FFTReal <float> &fft_second) const;

  ///
  /// Computes the crossing rate of the signal
  ///
    int compute_zcr(const std::vector<float> &x) const;


	///
	/// Returns the pitch (in Hz) of input frame x
	///
    float compute_pitch(std::vector<float> & x, FFTReal <float> &fft_first, FFTReal <float> &fft_second, unsigned int frame) ;
	
	///
	/// Returns true is the frame is unvoiced
	///
    bool unvoiced(float pot, float r1norm, float rmaxnorm, unsigned int *min, float zeros, float c0, unsigned int frame) ;
  
  ///
  /// Sets minimum
  ///
    void set_min(unsigned int min);

  public:
    PitchAnalyzer(	unsigned int fLen,			///< Frame length in samples
					unsigned int sFreq,			///< Sampling rate in Hertzs
					Window w=PitchAnalyzer::HAMMING,	///< Window type
					float min_F0 = MIN_F0,		///< Pitch range should be restricted to be above this value
					float max_F0 = MAX_F0		///< Pitch range should be restricted to be below this value
				 )
	  {
      frameLen = fLen;
      samplingFreq = sFreq;
      set_f0_range(min_F0, max_F0);
      set_window(w);
      pitch.resize(3);
    }

	///
    /// Operator (): computes the pitch for the given vector x
	///
    float operator()(const std::vector<float> & _x, FFTReal <float> &fft_first, FFTReal <float> &fft_second, unsigned int frame) {
      if (_x.size() != frameLen)
        return -1.0F;

      std::vector<float> x(_x); //local copy of input frame
      return compute_pitch(x, fft_first, fft_second, frame);
    }

	///
    /// Operator (): computes the pitch for the given "C" vector (float *).
    /// N is the size of the vector pointer by pt.
	///
    float operator()(const float * pt, unsigned int N, FFTReal <float> &fft_first, FFTReal <float> &fft_second, unsigned int frame) {
      if (N != frameLen)
        return -1.0F;

      std::vector<float> x(N); //local copy of input frame, size N
      std::copy(pt, pt+N, x.begin()); ///copy input values into local vector x
      return compute_pitch(x, fft_first, fft_second, frame);
    }

	///
    /// Operator (): computes the pitch for the given vector, expressed by the begin and end iterators
	///
    float operator()(std::vector<float>::const_iterator begin, std::vector<float>::const_iterator end, FFTReal <float> &fft_first, FFTReal <float> &fft_second, unsigned int frame) {

      if (end-begin != frameLen)
        return -1.0F;

      std::vector<float> x(end-begin); //local copy of input frame, size N
      std::copy(begin, end, x.begin()); //copy input values into local vector x
      return compute_pitch(x, fft_first,fft_second, frame);
    }
    
	///
    /// Sets pitch range: takes min_F0 and max_F0 in Hz, sets npitch_min and npitch_max in samples
	///
    void set_f0_range(float min_F0, float max_F0);
  };
}
#endif
