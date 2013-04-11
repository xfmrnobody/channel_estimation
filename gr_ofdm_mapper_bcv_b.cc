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

#include <gr_ofdm_mapper_bcv_b.h>
#include <gr_io_signature.h>
#include <stdexcept>
#include <string.h>
static int out_count=0;

gr_ofdm_mapper_bcv_b_sptr
gr_make_ofdm_mapper_bcv_b (const std::vector<gr_complex> &constellation, unsigned int msgq_limit, 
			 unsigned int occupied_carriers,const std::vector<gr_complex> &known_symbol, unsigned int fft_length)
{
  return gr_ofdm_mapper_bcv_b_sptr (new gr_ofdm_mapper_bcv_b (constellation, msgq_limit, 
							  occupied_carriers, known_symbol,fft_length));
}

// Consumes 1 packet and produces as many OFDM symbols of fft_length to hold the full packet
gr_ofdm_mapper_bcv_b::gr_ofdm_mapper_bcv_b (const std::vector<gr_complex> &constellation, unsigned int msgq_limit, 
					unsigned int occupied_carriers,const std::vector<gr_complex> &known_symbol, unsigned int fft_length)
  : gr_sync_block ("ofdm_mapper_bcv_b",
		   gr_make_io_signature (0, 0, 0),
		   gr_make_io_signature2 (1, 2, sizeof(gr_complex)*fft_length, sizeof(char))),
    d_constellation(constellation),
    d_msgq(gr_make_msg_queue(msgq_limit)), d_msg_offset(0), d_eof(false),
    d_occupied_carriers(occupied_carriers),
	d_known_symbol(known_symbol),
    d_fft_length(fft_length),
    d_bit_offset(0),
    d_pending_flag(0),
	d_pilot_flag(0),
    d_resid(0),
    d_nresid(0)
{
  if (!(d_occupied_carriers <= d_fft_length))
    throw std::invalid_argument("gr_ofdm_mapper_bcv_b: occupied carriers must be <= fft_length");

  // this is not the final form of this solution since we still use the occupied_tones concept,
  // which would get us into trouble if the number of carriers we seek is greater than the occupied carriers.
  // Eventually, we will get rid of the occupied_carriers concept.
  std::string carriers = "FE7F";

  // A bit hacky to fill out carriers to occupied_carriers length
  int diff = (d_occupied_carriers - 4*carriers.length()); 
  while(diff > 7) {
    carriers.insert(0, "f");
    carriers.insert(carriers.length(), "f");
    diff -= 8;
  }                            //使正中间两个子载波为0,即"...ffffFE7Fffff..."

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

    diff_left = (int)ceil((float)diff/2.0f);   // number of carriers to put on the left side
    c[0] = abc[(1 << diff_left) - 1];          // convert to bits and move to ASCI integer
    carriers.insert(0, c);
    
    diff_right = diff - diff_left;	       // number of carriers to put on the right side
    c[0] = abc[0xF^((1 << diff_right) - 1)];   // convert to bits and move to ASCI integer
        carriers.insert(carriers.length(), c);
  }  //作用与上面类似,针对occupide_length非8倍数,将上一步处理后剩余的子载波位置进行补1  最好是8的倍数
  
  // find out how many zeros to pad on the sides; the difference between the fft length and the subcarrier
  // mapping size in chunks of four. This is the number to pack on the left and this number plus any 
  // residual nulls (if odd) will be packed on the right. 
  diff = (d_fft_length/4 - carriers.length())/2; 
  //确定第几个子载波可以用，把它的序号写入d_subcarrier_map
  unsigned int i,j,k;
  for(i = 0; i < carriers.length(); i++) {
    char c = carriers[i];                            // get the current hex character from the string
    for(j = 0; j < 4; j++) {                         // walk through all four bits
      k = (strtol(&c, NULL, 16) >> (3-j)) & 0x1;     // convert to int and extract next bit
      if(k) {                                        // if bit is a 1, 
	d_subcarrier_map.push_back(4*(i+diff) + j);  // use this subcarrier
      }
    }
  }

  // make sure we stay in the limit currently imposed by the occupied_carriers
  if(d_subcarrier_map.size() > d_occupied_carriers) {
    throw std::invalid_argument("gr_ofdm_mapper_bcv_b: subcarriers allocated exceeds size of occupied carriers");
  }
  
  d_nbits = (unsigned long)ceil(log10(d_constellation.size()) / log10(2.0));  //对bpsk来说,为1
}

gr_ofdm_mapper_bcv_b::~gr_ofdm_mapper_bcv_b(void)
{
}

int gr_ofdm_mapper_bcv_b::randsym()
{
  return (rand() % d_constellation.size());
}

int
gr_ofdm_mapper_bcv_b::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
{
  gr_complex *out = (gr_complex *)output_items[0];
  
  unsigned int i=0;

  //printf("OFDM BPSK Mapper:  ninput_items: %d   noutput_items: %d\n", ninput_items[0], noutput_items);

  if(d_eof) {
    return -1;
  }
  
  if(!d_msg) /*判断是否是第一个包*/{
    d_msg = d_msgq->delete_head();	   // block, waiting for a message
    d_msg_offset = 0;
    d_bit_offset = 0;
    d_pending_flag = 1;			   // new packet, write start of packet flag
	d_pilot_flag=0;
    
    if((d_msg->length() == 0) && (d_msg->type() == 1)) {
      d_msg.reset();
      return -1;		// We're done; no more messages coming.
    }
  }

  char *out_flag = 0;
  if(output_items.size() == 2)
    out_flag = (char *) output_items[1];
  

  // Build a single symbol:
  // Initialize all bins to 0 to set unused carriers
  memset(out, 0, d_fft_length*sizeof(gr_complex));
  
  i = 0;

  unsigned char bits= 0;
  //以下,每8个bit发送一次,对应一个message,对子载波为8的倍数的bpsk来言,没有多余bit
  //while((d_msg_offset < d_msg->length()) && (i < d_occupied_carriers)) {
  while((d_msg_offset < d_msg->length()) && (i < d_subcarrier_map.size()/*除去occupied中间的两个子载波*/)) {
	
	//插入导频
	  if(d_pilot_flag%10==0){
		  for(i;i<d_subcarrier_map.size();i++){
		  out[d_subcarrier_map[i]]=d_known_symbol[i];
		  out_count++;
		  }
		  //printf("\n");
		  //printf("insert 1 pilot\n");
		  d_pilot_flag++;
		  continue;
	  }
		  
    // need new data to process
    if(d_bit_offset == 0) {
      d_msgbytes = d_msg->msg()[d_msg_offset];
      //printf("mod message byte: %x\n", d_msgbytes);
    }
    
    if(d_nresid > 0)/*对子载波数8的倍数bpsk来说,这一部分(处理上一个message剩余的比特)没有用*/ {
      // take the residual bits, fill out nbits with info from the new byte, and put them in the symbol
      d_resid |= (((1 << d_nresid)-1) & d_msgbytes) << (d_nbits - d_nresid);
      bits = d_resid;

      out[d_subcarrier_map[i]] = d_constellation[bits];
      i++;

      d_bit_offset += d_nresid;
      d_nresid = 0;
      d_resid = 0;
      printf("mod bit(r): %x   resid: %x   nresid: %d    bit_offset: %d\n", 
           bits, d_resid, d_nresid, d_bit_offset);
    }
    else {
      if((8 - d_bit_offset) >= d_nbits) {  // test to make sure we can fit nbits
	// take the nbits number of bits at a time from the byte to add to the symbol
	bits = ((1 << d_nbits)-1) & (d_msgbytes >> d_bit_offset);
	d_bit_offset += d_nbits;
	
	out[d_subcarrier_map[i]] = d_constellation[bits];
	out_count++;
	i++;
	if(i==d_subcarrier_map.size()){
		d_pilot_flag++;
		//printf("carrier NO:%d  \n",d_pilot_flag);
	}
      }
      else {  // if we can't fit nbits, store them for the next 
	// saves d_nresid bits of this message where d_nresid < d_nbits
		  //计算这一message剩余比特
	unsigned int extra = 8-d_bit_offset;
	d_resid = ((1 << extra)-1) & (d_msgbytes >> d_bit_offset);
	d_bit_offset += extra;
	d_nresid = d_nbits - extra;
      }
      
    }
            
    if(d_bit_offset == 8) {
	  //printf("next");
      d_bit_offset = 0;
      d_msg_offset++;
    }
  }

  // Ran out of data to put in symbol
  if (d_msg_offset == d_msg->length()) {
	  //printf("over\n");
    if(d_nresid > 0) {
      d_resid |= 0x00;
      bits = d_resid;
      d_nresid = 0;
      d_resid = 0;
    }
     //若发送比特数小于占用子载波数,则将剩余的占用子载波补上随机数
    //while(i < d_occupied_carriers) {   // finish filling out the symbol
    while(i < d_subcarrier_map.size()) {   // finish filling out the symbol
      out[d_subcarrier_map[i]] = d_constellation[randsym()];
      i++;
    }

    if (d_msg->type() == 1)	        // type == 1 sets EOF
      d_eof = true;
    d_msg.reset();   			// finished packet, free message
    assert(d_bit_offset == 0);
  }

  if (out_flag)
    out_flag[0] = d_pending_flag;
  d_pending_flag = 0;
  //printf("out_count:%d\n",out_count);

  return 1;  // produced symbol
}
