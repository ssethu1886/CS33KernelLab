#ifndef STENCIL_H
#define STENCIL_H

#include <stdint.h>
#include <stdio.h>

//Define all the test cases
#define T1I  10000
#define T1J  10000
#define T1K      1
#define T1S      0

#define T2I  128
#define T2J  128
#define T2K  128
#define T2S    8

#define T3I   64
#define T3J   64
#define T3K   64
#define T3S   20

#define T4I       4
#define T4J       4
#define T4K 1048576
#define T4S       2


#define NUM_TESTS (4)

void move_one_value(int i1, int j1, int i2, int j2, int Ni, int Nj, char In[Ni][Nj], char Out[Ni][Nj]);

void set_to_zero(int i, int j, int k, int Ni, int Nj, int Nk,float array[Ni][Nj][Nk]);

void macc_element(const float* In, float* Out, const float* stencil);

void reticulate_splines(do_all);

void compute_transpose(int Ni, int Nj,
                    const char In[Ni][Nj], char Out[Nj][Ni]);

void transpose_check(int Ni, int Nj,
                    const char In[Ni][Nj], char Out[Nj][Ni]);

void stencil_check(int Ni, int Nj, int Nk, int S, 
            const float In[Ni+S][Nj+S][Nk+S], float Out[Ni][Nj][Nk], 
            const float Stencil[S][S][S]);

void compute_stencil(int Ni, int Nj, int Nk, int S, 
            const float In[Ni+S][Nj+S][Nk+S], float Out[Ni][Nj][Nk], 
            const float Stencil[S][S][S]);

#endif