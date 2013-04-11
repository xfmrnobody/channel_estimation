/* -*- c++ -*- */
/*
 * Copyright 2006,2007,2008 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gr_ofdm_frame_acquisition_b1.h>
#include <gr_io_signature.h>
#include <gr_expj.h>
#include <gr_math.h>
#include <cstdio>
#include <fstream>

#define VERBOSE 0
#define M_TWOPI (2*M_PI)
#define MAX_NUM_SYMBOLS 1000

gr_ofdm_frame_acquisition_b1_sptr
gr_make_ofdm_frame_acquisition_b1 (unsigned int occupied_carriers, unsigned int fft_length, 
				unsigned int cplen,
				const std::vector<gr_complex> &known_symbol,
				const std::vector<gr_complex> &known_pilot,
				unsigned int max_fft_shift_len)
{
  return gr_ofdm_frame_acquisition_b1_sptr (new gr_ofdm_frame_acquisition_b1 (occupied_carriers, fft_length, cplen,
									known_symbol, known_pilot,max_fft_shift_len));
}

gr_ofdm_frame_acquisition_b1::gr_ofdm_frame_acquisition_b1 (unsigned occupied_carriers, unsigned int fft_length, 
						      unsigned int cplen,
						      const std::vector<gr_complex> &known_symbol,
						      const std::vector<gr_complex> &known_pilot,
						      unsigned int max_fft_shift_len)
  : gr_block ("ofdm_frame_acquisition",
	      gr_make_io_signature2 (2, 2, sizeof(gr_complex)*fft_length, sizeof(char)*fft_length),
	      gr_make_io_signature2 (2, 2, sizeof(gr_complex)*occupied_carriers, sizeof(char))),
    d_occupied_carriers(occupied_carriers),
    d_fft_length(fft_length),
    d_cplen(cplen),
    d_freq_shift_len(max_fft_shift_len),
    d_known_symbol(known_symbol),
	d_known_pilot(known_pilot),
    d_coarse_freq(0),
    d_phase_count(0),
	d_pilot_flag(0)
{
  d_symbol_phase_diff.resize(d_fft_length);
  d_known_phase_diff.resize(d_occupied_carriers);
  d_hestimate.resize(d_occupied_carriers);
  d_pilot.resize(d_occupied_carriers);
  /*******************************得到subcarrier_map***************************************/
  std::string carriers = "FE7F";
	  
  // A bit hacky to fill out carriers to occupied_carriers length
  int diff = (d_occupied_carriers - 4*carriers.length()); 
  while(diff > 7) {
	  carriers.insert(0, "f");
	  carriers.insert(carriers.length(), "f");
	  diff -= 8;
  }
  
  // if there's extras left to be processed
  // divide remaining to put on either side of current map
  // all of this is done to stick with the concept of a carrier map string that
  // can be later passed by the user, even though it'd be cleaner to just do this
  // on the carrier map itself
  int diff_left=0;
  int diff_right=0;
	  
  // dictionary to convert from integers to ascii hex representation
  char abc[16] = {'0', '1', '2', '3', '4', '5', '6', '7', 
				  '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
  if(diff > 0) {
	  char c[2] = {0,0};
	  
	  diff_left = (int)ceil((float)diff/2.0f);  // number of carriers to put on the left side
	  c[0] = abc[(1 << diff_left) - 1];         // convert to bits and move to ASCI integer
	  carriers.insert(0, c);
	  
	  diff_right = diff - diff_left;	      // number of carriers to put on the right side
	  c[0] = abc[0xF^((1 << diff_right) - 1)];  // convert to bits and move to ASCI integer
	  carriers.insert(carriers.length(), c);
  }
	  
  // It seemed like such a good idea at the time...
  // because we are only dealing with the occupied_carriers
  // at this point, the diff_left in the following compensates
  // for any offset from the 0th carrier introduced
  unsigned int x,y,z;
  for(x = 0; x < (d_occupied_carriers/4)+diff_left; x++) {
	  char c = carriers[x];
	  for(y = 0; y < 4; y++) {
		  z = (strtol(&c, NULL, 16) >> (3-y)) & 0x1;
		  if(z) {
			  d_subcarrier_map.push_back(4*x + y - diff_left);
		  }
	  }
  }
/***********************************************************************/

  unsigned int i = 0, j = 0;

  std::fill(d_known_phase_diff.begin(), d_known_phase_diff.end(), 0);
  for(i = 0; i < d_known_symbol.size()-2; i+=2) {
    d_known_phase_diff[i] = norm(d_known_symbol[i] - d_known_symbol[i+2]);
  }
  
  d_phase_lut = new gr_complex[(2*d_freq_shift_len+1) * MAX_NUM_SYMBOLS];
  for(i = 0; i <= 2*d_freq_shift_len; i++) {
    for(j = 0; j < MAX_NUM_SYMBOLS; j++) {
      d_phase_lut[j + i*MAX_NUM_SYMBOLS] =  gr_expj(-M_TWOPI*d_cplen/d_fft_length*(i-d_freq_shift_len)*j);
    }
  }
}

gr_ofdm_frame_acquisition_b1::~gr_ofdm_frame_acquisition_b1(void)
{
  delete [] d_phase_lut;
}

void
gr_ofdm_frame_acquisition_b1::forecast (int noutput_items, gr_vector_int &ninput_items_required)
{
  unsigned ninputs = ninput_items_required.size ();
  for (unsigned i = 0; i < ninputs; i++)
    ninput_items_required[i] = 1;
}

gr_complex
gr_ofdm_frame_acquisition_b1::coarse_freq_comp(int freq_delta, int symbol_count)
{
  //  return gr_complex(cos(-M_TWOPI*freq_delta*d_cplen/d_fft_length*symbol_count),
  //	    sin(-M_TWOPI*freq_delta*d_cplen/d_fft_length*symbol_count));

  return gr_expj(-M_TWOPI*freq_delta*d_cplen/d_fft_length*symbol_count);

  //return d_phase_lut[MAX_NUM_SYMBOLS * (d_freq_shift_len + freq_delta) + symbol_count];
}

void
gr_ofdm_frame_acquisition_b1::correlate(const gr_complex *symbol, int zeros_on_left)
{
  unsigned int i,j;
  
  std::fill(d_symbol_phase_diff.begin(), d_symbol_phase_diff.end(), 0);
  for(i = 0; i < d_fft_length-2; i++) {
    d_symbol_phase_diff[i] = norm(symbol[i] - symbol[i+2]);
  }

  // sweep through all possible/allowed frequency offsets and select the best
  int index = 0;
  float max = 0, sum=0;
  for(i =  zeros_on_left - d_freq_shift_len; i < zeros_on_left + d_freq_shift_len; i++) {
    sum = 0;
    for(j = 0; j < d_occupied_carriers; j++) {
      sum += (d_known_phase_diff[j] * d_symbol_phase_diff[i+j]);
    }
    if(sum > max) {
      max = sum;
      index = i;
    }
  }
  
  // set the coarse frequency offset relative to the edge of the occupied tones
  d_coarse_freq = index - zeros_on_left;
}

void
gr_ofdm_frame_acquisition_b1::calculate_equalizer(const gr_complex *symbol, int zeros_on_left)
{
  unsigned int i=0;

  // Set first tap of equalizer
  d_hestimate[0] = d_known_symbol[0] / 
    (coarse_freq_comp(d_coarse_freq,1)*symbol[zeros_on_left+d_coarse_freq]);

  // set every even tap based on known symbol
  // linearly interpolate between set carriers to set zero-filled carriers
  // FIXME: is this the best way to set this?
  for(i = 2; i < d_occupied_carriers; i+=2) {
    d_hestimate[i] = d_known_symbol[i] / 
      (coarse_freq_comp(d_coarse_freq,1)*(symbol[i+zeros_on_left+d_coarse_freq]));
    d_hestimate[i-1] = (d_hestimate[i] + d_hestimate[i-2]) / gr_complex(2.0, 0.0);    
  }

  // with even number of carriers; last equalizer tap is wrong
  if(!(d_occupied_carriers & 1)) {
    d_hestimate[d_occupied_carriers-1] = d_hestimate[d_occupied_carriers-2];
  }

  if(VERBOSE) {
    printf( "Equalizer setting:\n");
    for(i = 0; i < d_occupied_carriers; i++) {
      gr_complex sym = coarse_freq_comp(d_coarse_freq,1)*symbol[i+zeros_on_left+d_coarse_freq];
      gr_complex output = sym * d_hestimate[i];
      printf("sym: %+.4f + j%+.4f  ks: %+.4f + j%+.4f  eq: %+.4f + j%+.4f  ==>  %+.4f + j%+.4f\n", 
	      sym .real(), sym.imag(),
	      d_known_symbol[i].real(), d_known_symbol[i].imag(),
	      d_hestimate[i].real(), d_hestimate[i].imag(),
	      output.real(), output.imag());
    }
    printf( "\n");
  }
}
void
gr_ofdm_frame_acquisition_b1::calculate_estimation(const gr_complex *symbol,
												unsigned int d_phase_count,int zeros_on_left)
{
	std::ofstream f;
	char *path="./Estimation_Liner_Block";
	f.open(path,std::ios::app);
	unsigned int i = 0;
	for(i=0; i < d_subcarrier_map.size(); i++) {
		d_pilot[d_subcarrier_map[i]] = d_known_pilot[i] / 
			(coarse_freq_comp(d_coarse_freq,d_phase_count)*(symbol[d_subcarrier_map[i]+zeros_on_left+d_coarse_freq]));   
	}
	for(i=0;i<d_subcarrier_map.size();i++){
		f<<d_pilot[i]<<".\n";
	}
	f.close();
	if(VERBOSE) {
		printf( "Estimation setting:\n");
		for(i = 0; i < d_occupied_carriers; i++) {
			gr_complex sym = coarse_freq_comp(d_coarse_freq,d_phase_count)*symbol[i+zeros_on_left+d_coarse_freq];
			gr_complex output = sym * d_hestimate[i];
			printf( "sym: %+.4f + j%+.4f  ks: %+.4f + j%+.4f  eq: %+.4f + j%+.4f  ==>  %+.4f + j%+.4f\n", 
					sym .real(), sym.imag(),
					d_known_pilot[i].real(), d_known_pilot[i].imag(),
					d_pilot[i].real(), d_pilot[i].imag(),
					output.real(), output.imag());
		}
		printf("\n");
	}
}

int
gr_ofdm_frame_acquisition_b1::general_work(int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
{
  const gr_complex *symbol = (const gr_complex *)input_items[0];
  const char *signal_in = (const char *)input_items[1];

  gr_complex *out = (gr_complex *) output_items[0];
  char *signal_out = (char *) output_items[1];
  
  int unoccupied_carriers = d_fft_length - d_occupied_carriers;
  int zeros_on_left = (int)ceil(unoccupied_carriers/2.0);

  if(signal_in[0]) {
    d_phase_count = 1;
	d_pilot_flag=0;
    correlate(symbol, zeros_on_left);
    calculate_equalizer(symbol, zeros_on_left); 
    signal_out[0] = 1;
	for(unsigned int i = 0; i < d_occupied_carriers; i++) {
		out[i] = d_hestimate[i]*coarse_freq_comp(d_coarse_freq,d_phase_count)
				 *symbol[i+zeros_on_left+d_coarse_freq];
	}
	
  }
  else {
    signal_out[0] = 0;
	if(d_pilot_flag%10==0){
		calculate_estimation(symbol,d_phase_count,zeros_on_left);
	}
	for(unsigned int i = 0; i < d_occupied_carriers; i++) {
		out[i] = d_pilot[i]*coarse_freq_comp(d_coarse_freq,d_phase_count)
				 *symbol[i+zeros_on_left+d_coarse_freq];
	}
	
	//printf("carrier NO:%d\n",d_pilot_flag);
	d_pilot_flag++;
	
  } 
  
  d_phase_count++;
  if(d_phase_count == MAX_NUM_SYMBOLS) {
    d_phase_count = 1;
  }

  consume_each(1);
  return 1;
}
