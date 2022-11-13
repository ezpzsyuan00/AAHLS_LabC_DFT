#include<math.h>
#include "dft.h"
#include "coefficients1024.h"
#include<stdio.h>
#include <stdlib.h>
#include<iostream>
/*
void dft1(DTYPE real_sample[1024], DTYPE imag_sample[1024],DTYPE real_op[1024],DTYPE imag_op[1024])
{
	//Write your code here
	#pragma HLS INTERFACE mode=s_axilite port=return
	#pragma HLS INTERFACE mode=m_axi bundle=A depth=1024 num_read_outstanding = 1 max_write_burst_length = 2 max_read_burst_length = 2 port=real_sample
	#pragma HLS INTERFACE mode=m_axi bundle=B depth=1024 num_read_outstanding = 1 max_write_burst_length = 2 max_read_burst_length = 2 port=imag_sample
	#pragma HLS INTERFACE mode=m_axi bundle=C depth=1024 num_write_outstanding = 2 max_write_burst_length = 2 max_read_burst_length = 2 port=real_op
	#pragma HLS INTERFACE mode=m_axi bundle=D depth=1024 num_write_outstanding = 2 max_write_burst_length = 2 max_read_burst_length = 2 port=imag_op


	dft_label0:for(int m=0; m<512; m++){
		int index1;
		int index2;
		DTYPE c;
		DTYPE s;
		DTYPE wc;
		DTYPE ws;
		DTYPE rco;
		DTYPE iso;
		DTYPE rso;
		DTYPE ico;
		DTYPE rce;
		DTYPE ise;
		DTYPE rse;
		DTYPE ice;
		DTYPE Ro;
		DTYPE Io;
		DTYPE Re;
		DTYPE Ie;
		dft_label1:for(int k=0; k<512; k++){
			#pragma HLS PIPELINE rewind
			index1 = (2*m*k)&1023;
			index2 = k;
			c = cos_coefficients_table[index1];
			s = sin_coefficients_table[index1];
			wc = cos_coefficients_table[index2];
			ws = sin_coefficients_table[index2];
			rce = real_sample[2*m]*c;
			ise = imag_sample[2*m]*s;
			rse = real_sample[2*m]*s;
			ice = imag_sample[2*m]*c;

			rco = real_sample[2*m+1]*c;
			iso = imag_sample[2*m+1]*s;
			rso = real_sample[2*m+1]*s;
			ico = imag_sample[2*m+1]*c;


			Ro = (rco+iso);
			Io = (rso-ico);
			Re = (rce+ise);
			Ie = (rse-ice);

			real_op[k] += (Re + (Ro*wc-Io*ws));
			imag_op[k] += (Ie + (Ro*ws+Io*wc));
			real_op[k+512] += (Re - (Ro*wc-Io*ws));
			imag_op[k+512] += (Ie - (Ro*ws+Io*wc));
		}
	}
}
*/

void dft(DTYPE real_sample[1024], DTYPE imag_sample[1024],DTYPE real_op[1024],DTYPE imag_op[1024])
{
	//Write your code here
	#pragma HLS INTERFACE mode=s_axilite port=return
	#pragma HLS INTERFACE mode=m_axi bundle=A depth=1024 port=real_sample
	#pragma HLS INTERFACE mode=m_axi bundle=B depth=1024 port=imag_sample
	#pragma HLS INTERFACE mode=m_axi bundle=C depth=1024 port=real_op
	#pragma HLS INTERFACE mode=m_axi bundle=D depth=1024 port=imag_op

	dft_label0:for(int i=0; i<1024; i++){
		int index;
		DTYPE c;
		DTYPE s;
		DTYPE rc;
		DTYPE is;
		DTYPE rs;
		DTYPE ic;

		dft_label1:for(int j=0; j<1024; j++){
			#pragma HLS PIPELINE II=1 rewind
			index = (i*j)&1023;
			c = cos_coefficients_table[index];
			s = sin_coefficients_table[index];
			rc = real_sample[i]*c;
			is = imag_sample[i]*s;
			rs = real_sample[i]*s;
			ic = imag_sample[i]*c;

			real_op[j] += (rc+is);
			imag_op[j] += (rs-ic);
		}
	}
}

