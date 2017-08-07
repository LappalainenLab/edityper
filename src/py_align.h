/*
 * py_align.h for program py_align
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <utility>

#define DEBUG 0

using namespace std;

string nw_align(string, string, int);

string nw_free_trail(string, string);

string nw_align_mem(string, string, int, int, int);

string nw_align_aff(string, string, int, int);

string nw_align_aff_mem(string, string, int, int, int, int);

string nw_align_aff_param(string, string, int, int, int, int);

int max_2(int, int);
int max_3(int, int, int);

void  dpm_init( int **, int, int);
void  P_init( int **, int, int);
void  Q_init( int **, int, int);
