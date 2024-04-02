#ifndef HEADERS_H
#define HEADERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "gnuplot_i.h"
#include "matrix2.h"
#include "image.h"

#define DATA_DIR  "data/"
#define ADATA_LEN 60000.0
#define BDATA_LEN 140000.0

typedef struct gauss_t {
    int id;
    VEC *mu;
    MAT *sigma;
    MAT *dataset;
    double prior;
} gauss_t;

typedef double (*Disc)(VEC *, gauss_t *);

typedef struct batch_t {
    char bname[32];
    size_t n;
    size_t c;
    size_t d;
    
    Disc g;    
    FILE *fdata;
    gauss_t **dist;
} batch_t;

#endif // HEADERS_H
