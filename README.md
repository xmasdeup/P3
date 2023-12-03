PAV - P3: estimación de pitch
=============================

Esta práctica se distribuye a través del repositorio GitHub [Práctica 3](https://github.com/albino-pav/P3).
Siga las instrucciones de la [Práctica 2](https://github.com/albino-pav/P2) para realizar un `fork` de la
misma y distribuir copias locales (*clones*) del mismo a los distintos integrantes del grupo de prácticas.

Recuerde realizar el *pull request* al repositorio original una vez completada la práctica.

Ejercicios básicos
------------------

- Complete el código de los ficheros necesarios para realizar la estimación de pitch usando el programa
  `get_pitch`.

   * Complete el cálculo de la autocorrelación e inserte a continuación el código correspondiente.

          void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    
          for (unsigned int l = 0; l < r.size(); ++l) 
            
            {
            r[l] = 0;
            for (unsigned int m = 0; m < (x.size()-l); ++m){
              r[l] += x[m]*x[m+l];
              }
              r[l] = r[l]/x.size();
            }
          
            if (r[0] == 0.0F) 

              r[0] = 1e-10; 
            } 

   * Inserte una gŕafica donde, en un *subplot*, se vea con claridad la señal temporal de un segmento de
     unos 30 ms de un fonema sonoro y su periodo de pitch; y, en otro *subplot*, se vea con claridad la
	 autocorrelación de la señal y la posición del primer máximo secundario.

	 NOTA: es más que probable que tenga que usar Python, Octave/MATLAB u otro programa semejante para
	 hacerlo. Se valorará la utilización de la biblioteca matplotlib de Python.

   * Determine el mejor candidato para el periodo de pitch localizando el primer máximo secundario de la
     autocorrelación. Inserte a continuación el código correspondiente.
    
          unsigned int PitchAnalyzer::amdf(const vector<float> &x, vector<float> &distance) const{
            
            float_t minIndex = 1000;

            float_t mean = 0;

            unsigned int index = 0;

            for (unsigned int lag = npitch_min-2; lag <= distance.size() ;lag++){

                distance[lag] = 0;

              for(unsigned int n = 0; n < (x.size() -lag); ++n)
              {

                distance[lag] += abs(x[n]-x[n+lag]);

              }

              distance[lag] = distance[lag]/(x.size()-lag);

              if(lag<x.size()/2) mean+=distance[lag];
              
              if(lag<x.size()/2) ;

                if(minIndex>distance[lag]) 
                {

                  minIndex = distance[lag];
                  index = lag;

                }

            }

            if(minIndex < 1.5e-4) index = 0; 

          return index;
         }

          unsigned int min = amdf(x,distance);

   * Implemente la regla de decisión sonoro o sordo e inserte el código correspondiente.

   * Puede serle útil seguir las instrucciones contenidas en el documento adjunto `código.pdf`.

- Una vez completados los puntos anteriores, dispondrá de una primera versión del estimador de pitch. El 
  resto del trabajo consiste, básicamente, en obtener las mejores prestaciones posibles con él.

  * Utilice el programa `wavesurfer` para analizar las condiciones apropiadas para determinar si un
    segmento es sonoro o sordo. 
	
	  - Inserte una gráfica con la estimación de pitch incorporada a `wavesurfer` y, junto a ella, los 
	    principales candidatos para determinar la sonoridad de la voz: el nivel de potencia de la señal
		(r[0]), la autocorrelación normalizada de uno (r1norm = r[1] / r[0]) y el valor de la
		autocorrelación en su máximo secundario (rmaxnorm = r[lag] / r[0]).

		Puede considerar, también, la conveniencia de usar la tasa de cruces por cero.

	    Recuerde configurar los paneles de datos para que el desplazamiento de ventana sea el adecuado, que
		en esta práctica es de 15 ms.

      - Use el estimador de pitch implementado en el programa `wavesurfer` en una señal de prueba y compare
	    su resultado con el obtenido por la mejor versión de su propio sistema.  Inserte una gráfica
		ilustrativa del resultado de ambos estimadores.
     
		Aunque puede usar el propio Wavesurfer para obtener la representación, se valorará
	 	el uso de alternativas de mayor calidad (particularmente Python).
  
  * Optimice los parámetros de su sistema de estimación de pitch e inserte una tabla con las tasas de error
    y el *score* TOTAL proporcionados por `pitch_evaluate` en la evaluación de la base de datos 
	`pitch_db/train`..

Ejercicios de ampliación
------------------------

- Usando la librería `docopt_cpp`, modifique el fichero `get_pitch.cpp` para incorporar los parámetros del
  estimador a los argumentos de la línea de comandos.
  
  Esta técnica le resultará especialmente útil para optimizar los parámetros del estimador. Recuerde que
  una parte importante de la evaluación recaerá en el resultado obtenido en la estimación de pitch en la
  base de datos.

  * Inserte un *pantallazo* en el que se vea el mensaje de ayuda del programa y un ejemplo de utilización
    con los argumentos añadidos.

- Implemente las técnicas que considere oportunas para optimizar las prestaciones del sistema de estimación
  de pitch.

  Entre las posibles mejoras, puede escoger una o más de las siguientes:

  * Técnicas de preprocesado: filtrado paso bajo, diezmado, *center clipping*, etc.

    - Filtrado paso bajo:

          void ButterWorthFilter::applyFilter(float spectrum [], const float filter [])
          {

          for(int i = 0; i< samples; i++)
            {

            spectrum[i] *= filter[i];

              spectrum[i+samples] *= filter[i];

            }
          }

    - Center clipping:

                void ButterWorthFilter::center_clipping(std::vector<float> &signal, int rate) const
          {

            float_t interval = rate/10;

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

              threshold = maximum * 0.20;

              for(int32_t k = 0; k < interval; k++)
              {

                if( *(iX+k) >= (threshold)) ;

                else if(*(iX+k) <= (-threshold)) ; 
                else 
                {
                  *(iX+k) = 0;
                }
              }

            }

            if(iX < signal.end()){

              int missing = signal.end() - iX;

              for(int32_t k = 0; k < missing; k++)

              {
              if( *(iX+k) >= (threshold)) ;

              else if(*(iX+k) <= (-threshold)) ; 

              else 

              {

                *(iX+k) = 0;

              }

              }

            }

          }

  * Técnicas de postprocesado: filtro de mediana, *dynamic time warping*, etc.
  * Métodos alternativos a la autocorrelación: procesado cepstral, *average magnitude difference function*
    (AMDF), etc.
  
    - Procesado cepstral:

          void PitchAnalyzer::calcularCepstrum(const std::vector<float> &signal, std::vector<float> &cepstrum, int size_fft, FFTReal <float> &fft_first, FFTReal <float> &fft_second)const
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

            fft_second.rescale(out_cepstrum);

            for(int i = 0; i<(size_fft/4); i++)
            {
              cepstrum[i] = pow(out_cepstrum[i],2) + pow(out_cepstrum[i+size_fft/4],2);
              // cout<<cepstrum[i]<<"\n";
            }

          }

    - Average magnitude difference function (AMFD): 

          unsigned int PitchAnalyzer::amdf(const vector<float> &x, vector<float> &distance) const{
          
          float_t minIndex = 1000;

          unsigned int index = 0;

          std::cout<<MAX_F0<<"\n";

          for (unsigned int lag = npitch_min; lag <= distance.size() ;lag++)

          {

              distance[lag] = 0;

            for(unsigned int n = 0; n < (x.size() -lag); ++n)

            {

              distance[lag] += abs(x[n]-x[n+lag]);

            }

            distance[lag] = distance[lag]/(x.size()-lag);

            if(lag<x.size()/2) ;

            if(minIndex>distance[lag]) 
            {

              minIndex = distance[lag];
              index = lag;
              
            }

          }
          if(minIndex < 1.5e-4) index = 0; 

          std::cout<<minIndex<<"\t";

          std::cout<<index<<"\n";

          std::cout<<"End Of Frame\n";

          return index;
           }

  * Optimización **demostrable** de los parámetros que gobiernan el estimador, en concreto, de los que
    gobiernan la decisión sonoro/sordo.
  * Cualquier otra técnica que se le pueda ocurrir o encuentre en la literatura.

  Encontrará más información acerca de estas técnicas en las [Transparencias del Curso](https://atenea.upc.edu/pluginfile.php/2908770/mod_resource/content/3/2b_PS%20Techniques.pdf)
  y en [Spoken Language Processing](https://discovery.upc.edu/iii/encore/record/C__Rb1233593?lang=cat).
  También encontrará más información en los anexos del enunciado de esta práctica.

  Incluya, a continuación, una explicación de las técnicas incorporadas al estimador. Se valorará la
  inclusión de gráficas, tablas, código o cualquier otra cosa que ayude a comprender el trabajo realizado.

  También se valorará la realización de un estudio de los parámetros involucrados. Por ejemplo, si se opta
  por implementar el filtro de mediana, se valorará el análisis de los resultados obtenidos en función de
  la longitud del filtro.
   

Evaluación *ciega* del estimador
-------------------------------

Antes de realizar el *pull request* debe asegurarse de que su repositorio contiene los ficheros necesarios
para compilar los programas correctamente ejecutando `make release`.

Con los ejecutables construidos de esta manera, los profesores de la asignatura procederán a evaluar el
estimador con la parte de test de la base de datos (desconocida para los alumnos). Una parte importante de
la nota de la práctica recaerá en el resultado de esta evaluación.
