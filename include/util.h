#ifndef UTIL_H
#define UTIL_H

#include "headers.h"

double ranf(double m);
double rang(double mu, double sigma);
VEC *rang2D(VEC *xy, VEC *mu, MAT *sigma);

double m_det(MAT *m);
double gauss_eval(gauss_t *g, VEC *x);
double bhatta_err(gauss_t *c1, gauss_t *c2);
void sample_mean(MAT *data, VEC *out);
void sample_cov(MAT *data, VEC *mean, MAT *out);
void compute_mle(gauss_t *g, size_t n);

FILE *setup_dataset(char *fname, gauss_t *dist, size_t n);
MAT *gen_dataset(FILE *fp, MAT *m, gauss_t *dist, size_t n);
size_t load_dataset(FILE *fp, MAT *data);
size_t trim_zeros(MAT *m);
void norm_sort(MAT *m);
int dataset_len_cmp(const void *a, const void *b);

MAT *image_rgmat(Image *img, MAT *m);
MAT *image_ycbcrmat(Image *img, MAT *m);

#endif // UTIL_H