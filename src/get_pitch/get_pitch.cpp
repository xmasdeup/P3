/// @file

#include <iostream>
#include <fstream>
#include <string.h>
#include <errno.h>

#include "wavfile_mono.h"
#include "pitch_analyzer.h"
#include "FFTReal.h"
#include	"def.h"
#include	"DynArray.h"
#include	"OscSinCos.h"

#include "docopt.h"

#define FRAME_LEN   0.030 /* 30 ms. */
#define FRAME_SHIFT 0.015 /* 15 ms. */

using namespace std;
using namespace upc;
using namespace ffft;

static const char USAGE[] = R"(
get_pitch - Pitch Estimator 

Usage:
    get_pitch [options] <input-wav> <output-txt>
    get_pitch (-h | --help)
    get_pitch --version

Options:
    -h, --help  Show this screen
    --version   Show the version of the project

Arguments:
    input-wav   Wave file with the audio signal
    output-txt  Output file: ASCII file with the result of the estimation:
                    - One line per frame with the estimated f0
                    - If considered unvoiced, f0 must be set to f0 = 0
)";

int main(int argc, const char *argv[]) {
	/// \TODO 
	///  Modify the program syntax and the call to **docopt()** in order to
	///  add options and arguments to the program.
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
        {argv + 1, argv + argc},	// array of arguments, without the program name
        true,    // show help if requested
        "2.0");  // version string

	std::string input_wav = args["<input-wav>"].asString();
	std::string output_txt = args["<output-txt>"].asString();

  // Read input sound file
  unsigned int rate;
  vector<float> x;
  if (readwav_mono(input_wav, rate, x) != 0) {
    cerr << "Error reading input file " << input_wav << " (" << strerror(errno) << ")\n";
    return -2;
  }
  // HACER NORMALIZACION DE X
  int n_len = rate * FRAME_LEN;
  int n_shift = rate * FRAME_SHIFT;
  int order = 4;
  int cutoffFrequency = 1100;
  int16_t DFTsamples = 1024;
  int16_t f_inc = round(rate/DFTsamples);
  int size_fft = 0;

  for(int size_f = 1; size_f < n_len; size_f *=2) size_fft = size_f;


  // Define analyzer
  PitchAnalyzer analyzer(n_len, rate, PitchAnalyzer::RECT, 100, 550);
  ButterWorthFilter butter(DFTsamples/2);
  FFTReal <float> fft_caller(DFTsamples);

  FFTReal <float> fft_first(size_fft);
  FFTReal <float> fft_second(size_fft/2);
  
  vector<float>::iterator iF;
  float filter[DFTsamples];
  float spectrum[DFTsamples];  

  for(int i = 1; i < DFTsamples/2; i++){

    float_t freq = i * f_inc;
    double denum = 1 + 0.31*pow(freq/cutoffFrequency,2*order);
    filter[i] = 1.0/sqrt(denum);

  }
    filter[0] = 1;
    // for(int i = 0; i<DFTsamples/2; i++)
    //  {
    //    cout<<i<<": ";
    //    cout<<filter[i]<<"\n";
    //  }

  for (iF = x.begin(); iF + DFTsamples < x.end(); iF = iF + DFTsamples) {
    
    float dft[DFTsamples]; //local copy of input frame, size N
    float sum = 0;

    for(uint16_t i = 0; i< DFTsamples; i++)
    {
      dft[i] = *(iF +i); 
      // cout<<i<<": ";
      // cout<<dft[i]<<"\n";
      sum += dft[i];
    }
    
    fft_caller.do_fft(spectrum,dft);

    butter.applyFilter(spectrum,filter);

    //  for(int i = 0; i<DFTsamples; i++)
    //  {
    //    cout<<i<<": ";
    //    cout<<spectrum[i]<<"\n";
    //  }

    fft_caller.do_ifft(spectrum,dft);

    fft_caller.rescale(dft);

    for(uint16_t i = 0; i<DFTsamples; i++)
    {
      *(iF+i) = dft[i];
      //cout<<(dft[i])<<"\n";

    }
    //cout<<"END OF FRAME\n\n\n\n\n";
  }
  
  butter.center_clipping(x,rate);


  /// \TODO
  /// Preprocess the input signal in order to ease pitch estimation. For instance,
  /// central-clipping or low pass filtering may be used.
  
  // Iterate for each frame and save values in f0 vector
  vector<float>::iterator iX;
  vector<float> f0;
  for (iX = x.begin(); iX + n_len < x.end(); iX = iX + n_shift) {
    float f = analyzer(iX, iX + n_len, fft_first, fft_second);
    f0.push_back(f);
  }

  /// \TODO
  /// Postprocess the estimation in order to supress errors. For instance, a median filter
  /// or time-warping may be used.

  // Write f0 contour into the output file
  ofstream os(output_txt);
  if (!os.good()) {
    cerr << "Error reading output file " << output_txt << " (" << strerror(errno) << ")\n";
    return -3;
  }

  os << 0 << '\n'; //pitch at t=0
  for (iX = f0.begin(); iX != f0.end(); ++iX) 
    os << *iX << '\n';
  os << 0 << '\n';//pitch at t=Dur

  return 0;
}
