#include "util.h"

double ranf(double m) {
    return (m * rand() / (double)RAND_MAX);
}

double rang(double mu, double sigma) {
    static double x1, x2, y1, y2;
    double w;

    do {
        x1 = 2.0 * ranf(1) - 1;
        x2 = 2.0 * ranf(1) - 1;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);

    y1 = x1 * w;
    y2 = x2 * w;

    return y1 * sigma + mu;
}

VEC *rang2D(VEC *xy, VEC *mu, MAT *sigma) {
    xy->ve[0] = rang(v_get_val(mu, 0), m_get_val(sigma, 0, 0));
    xy->ve[1] = rang(v_get_val(mu, 1), m_get_val(sigma, 1, 1));

    return xy;
}

double gauss_eval(gauss_t *g, VEC *x) {
    double a, b, c;
    VEC *y, *z;
    MAT *inv;

    inv = m_inverse(g->sigma, MNULL);
    y = v_sub(x, g->mu, VNULL);
    z = mv_mlt(inv, y, VNULL);

    a = 1 / (pow(2 * M_PI, (double)x->dim / 2) * sqrt(m_det(g->sigma)));
    b = -0.5 * in_prod(y, z);
    c = a * exp(b);

    v_free(y); v_free(z);
    m_free(inv);

    return c;
}

double bhatta_err(gauss_t *c1, gauss_t *c2) {
    MAT *s, *si;
    VEC *mu_diff;
    double ds, d1, d2, e, k;

    mu_diff = v_sub(c2->mu, c1->mu, VNULL);

    s = m_add(c1->sigma, c2->sigma, MNULL);
    sm_mlt(0.5, s, s);

    si = m_inverse(s, MNULL);

    ds = m_det(s);
    d1 = m_det(c1->sigma); d2 = m_det(c2->sigma);

    k = 0.125 * in_prod(mu_diff, mv_mlt(si, mu_diff, VNULL));
    k += 0.5 * log(ds / sqrt(d1 * d2));  

    e = pow(c1->prior, 0.5) * pow(c2->prior, 0.5) * exp(-k);

    m_free(s); m_free(si);
    v_free(mu_diff);

    return e;
}

void sample_mean(MAT *data, VEC *out) {
    MAT *tdata = m_transp(data, MNULL);
    VEC *avg = v_get(data->m);
    size_t i;

    for (i = 0; i < avg->dim; i += 1) avg->ve[i] = 1.0 / data->m;

    mv_mlt(tdata, avg, out);

    v_free(avg);
    m_free(tdata);
}

void sample_cov(MAT *data, VEC *mean, MAT *out) {
    MAT *tdata = m_get(data->n, data->m);
    size_t i, k;

    for (i = 0; i < data->n; i += 1) {
        for (k = 0; k < data->m; k += 1) {
            tdata->me[i][k] = data->me[k][i] - mean->ve[i];
        }
    }

    m_mlt(tdata, data, out);
    sm_mlt(1.0 / (data->m - 1), out, out);

    m_free(tdata);
}

void compute_mle(gauss_t *g, size_t n) {
    MAT *mle = m_copy(g->dataset, MNULL);

    mle = m_resize(mle, n, 2);

    sample_mean(mle, g->mu);
    sample_cov(mle, g->mu, g->sigma);

    m_free(mle);
}

FILE *setup_dataset(char *fname, gauss_t *dist, size_t n) {
    char path[MAX_FPATH];
    FILE *fp;

    sprintf(path, "%s%s", DATA_DIR, fname);

    if ((fp = fopen(path, "r+"))) {
        if (!load_dataset(fp, dist->dataset)) {
            printf("Failed to load data matrix '%s'...\n", path);
            fclose(fp);
        }
    } else if ((fp = fopen(path, "w+"))) {
        dist->dataset = gen_dataset(fp, dist->dataset, dist, n);
    } else {
        printf("Failed to open '%s'...\n", path);
    }

    return fp;
}

MAT *gen_dataset(FILE *fp, MAT *m, gauss_t *dist, size_t n) {
    VEC *v;
    MAT *out;

    if (!fp) return NULL;
    if (!m) out = m_resize(MNULL, n, 2);
    else out = m_resize(m, n, 2);

    v = v_get(2);

    while (n--) {
        rang2D(v, dist->mu, dist->sigma);
        set_row(out, n, v);
    }

    m_foutput(fp, out);

    v_free(v); 
    
    return out;
}

size_t load_dataset(FILE *fp, MAT *data) {
    fseek(fp, 0, SEEK_SET);
    
    data = m_finput(fp, data);

    return data->m;
}

size_t trim_zeros(MAT *m) {
    MAT *s;
    size_t l;
    VEC *v;
    v = v_get(m->n);

    norm_sort(m);

    for (l = m->m - 1; l >= 0; l -= 1) {
        get_row(m, l, v);
        if (v_norm2(v) == 0) break;
        set_row(m, m->m - 1 - l, v);
    }

    m_resize(m, m->m - 1 - l, m->n);

    return l;
}

void norm_sort(MAT *m) {
    MAT *s;
    VEC *v, *n;
    PERM *p;
    size_t i;

    p = px_get(m->m);
    n = v_get(m->m);
    v = v_get(m->n);

    for (i = 0; i < m->m; i += 1) {
        get_row(m, i, v);
        v_set_val(n, i, v_norm2(v));
    }

    v_sort(n, p);
    
    s = px_rows(p, m, MNULL);
    
    m_copy(s, m);

    m_free(s);
    v_free(v); v_free(n);
    px_free(p);
}

int dataset_len_cmp(const void *a, const void *b) {
    const MAT *l = (const MAT *)a;
    const MAT *r = (const MAT *)b;

    return (r->m - l->m);
}

double m_det(MAT *m) {
    MAT *a = m_copy(m, MNULL);
    PERM *p = px_get(a->m);
    double det;
    int i;

    det = 1;
    LUfactor(a, p);

    for (i = 0; i < a->m; i += 1) {
        det *= a->me[i][i];
    }

    for (i = 0; i < a->m; i += 1) {
        if (p->pe[i] != i) det = -det;
    }

    m_free(a);
    px_free(p);

    return det;
}

MAT *image_rgmat(Image *img, MAT *m) {
    MAT *out;
    size_t i, r;
    double n;

    if (!img) return NULL;
    if (!m) out = m_get(img->m * img->n, 2);
    else out = m_resize(m, img->m * img->n, 2);

    for (i = 0; i < img->size; i += 3) {
        for (r = 0, n = 0; r < 3; r += 1) n += img->data[i + r];
        if (n > 0) {
            m_set_val(out, i / 3, 0, (double)img->data[i] / n);
            m_set_val(out, i / 3, 1, (double)img->data[i + 1] / n);
        }
    }

    return out;
}