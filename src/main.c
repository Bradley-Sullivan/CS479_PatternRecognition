#include "headers.h"
#include "disc.h"
#include "util.h"

#define PLOT_DIR "plots/"

void mle_exp(MAT **datasets, batch_t *batch, int it);
size_t classify(batch_t *batch, MAT *data, int class);
void plot_batch(batch_t *b);
void plot_data(MAT **data, int n, char *fname, char *title);

void plot_roc(MAT *data, char *fname, char *title) {
    gnuplot_ctrl *plot;
    size_t i, k;
    char pbuf[32], *p;

    plot = gnuplot_init();

    gnuplot_cmd(plot, "set terminal pngcairo size 800,600 enhanced font 'Consolas,12'");
    gnuplot_cmd(plot, "set output '%s%s.png'", PLOT_DIR, fname);

    gnuplot_set_xlabel(plot, "FP");
    gnuplot_set_ylabel(plot, "FN");
    gnuplot_cmd(plot, "set title \"%s\"", title);

    gnuplot_cmd(plot, "set style line 1 lc rgb \"red\" linetype 1 linewidth 2 pointtype 7 pointsize 1.5 notitle");
    gnuplot_cmd(plot, "set style textbox opaque");

    gnuplot_cmd(plot, "plot '-' using 1:2 with linespoints linestyle 1,\\\n'' using 1:2:3 with labels center boxed notitle");

    for (i = 0; i < data->m; i += 1) {
        for (k = 0, p = pbuf; k < data->n; k += 1) {
            p += sprintf(p, "%lf ", data->me[i][k]);
        }
        gnuplot_cmd(plot, "%s", pbuf);
    }
    gnuplot_cmd(plot, "e\n");

    for (i = 0; i < data->m; i += 1) {
        for (k = 0, p = pbuf; k < data->n; k += 1) {
            p += sprintf(p, "%lf ", data->me[i][k]);
        }
        gnuplot_cmd(plot, "%s", pbuf);
    }
    gnuplot_cmd(plot, "e\n");

    gnuplot_close(plot);
}

void face_validate(Image *img, Image *ref, double *fpr, double *fnr) {
    size_t i, s, m, r;
    double fp, fn, tp, tn;

    fp = 0; fn = 0; tp = 0; tn = 0;
    for (i = 0, *fpr = 0, *fnr = 0; i < img->size; i += 3) {
        for (s = 0, m = 0, r = 0; s < 3; s += 1) {
            m += img->data[i + s];
            r += ref->data[i + s];
        }

        if (m && !r) {          // false positive
            fp += 1;
        } else if (!m && r) {   // false negative
            fn += 1;
        } else if (!m && !r) {  // true negative
            tn += 1;
        } else if (m && r) {    // true positive
            tp += 1;
        }
    }

    *fpr = fp / (fp + tn);
    *fnr = fn / (fn + tp);
}

void face_test(gauss_t *color, Image *out, Image *img, MAT *data, double thresh) {
    size_t k, s;
    double d;
    VEC *v;

    v = v_get(data->n);
    // iterate RG vectors
    for (k = 0; k < data->m; k += 1) {
        // get and test RG vectors against model
        get_row(data, k, v);
        d = gauss_eval(color, v);
        // if likelihood above threshold
        if (d >= thresh) {
            // copy over RGB values to output image
            for (s = 0; s < 3; s += 1) {
                out->data[k * 3 + s] = img->data[k * 3 + s];
            } 
        }
    }

    v_free(v);
}

void face_detect(gauss_t *color, char *ifname, char *rfname, size_t n, double step) {
    Image *img, *ref, *out;
    MAT *tdata, *roc;
    FILE *fp;
    VEC *v;
    double d, t;
    size_t i, k, s, e;
    double *fpr, *fnr, err;
    char fnbuf[MAX_FPATH], tbuf[MAX_FPATH];
    
    // allocate ROC stat buffers
    fpr = malloc(sizeof(size_t) * n);
    fnr = malloc(sizeof(size_t) * n);

    // load input images
    img = load_image(ifname);
    ref = load_image(rfname);
    out = new_image(img->m, img->n, img->q);

    printf("Getting '%s' RG vectors for testing...\n", ifname);
    // get input image dataset
    tdata = image_rgmat(img, MNULL);
    v = v_get(tdata->n);

    printf("Testing %llu batches with threshold-step of %lf...\n", n, step);
    // iterate batches
    for (i = 0, t = step, err = 999; i < n; i += 1, t += step) {
        printf("\tTesting batch %llu with threshold %lf...\n", i, t);
        face_test(color, out, img, tdata, t);

        // compute current threshold ROC stats
        face_validate(out, ref, fpr + i, fnr + i);
        if (fabs(fpr[i] - fnr[i]) <= err) { e = i; err = fabs(fpr[i] - fnr[i]); }

        // reset output image
        memset(out->data, 0, sizeof(uint8_t) * out->size);
    }

    // build ROC curve dataset -> [ fp, fn, threshold ]
    roc = m_get(n, 3);
    v_resize(v, 3);
    for (i = 0, t = 0; i < n; i += 1, t += step) {
        v_set_val(v, 0, fnr[i]);
        v_set_val(v, 1, fpr[i]);
        v_set_val(v, 2, t);
        set_row(roc, i, v);
    }

    face_test(color, out, img, tdata, m_get_val(roc, e, 2));
    sprintf(fnbuf, "err_thresh_%s", ifname);
    write_image(fnbuf, out);
    
    sprintf(fnbuf, "roc_data_%s.mat", ifname);
    fp = fopen(fnbuf, "w+");
    m_foutput(fp, roc);

    printf("Plotting ROC curve for '%s' test batch (n = %llu, step = %lf)...\n", ifname, n, step);
    // plot ROC curve
    sprintf(fnbuf, "roc_%s", ifname);
    sprintf(tbuf, "ROC for %s (n = %llu, step = %.02lf)\n", ifname, n, step);
    plot_roc(roc, fnbuf, tbuf);

    del_image(img); del_image(out); del_image(ref);
    m_free(tdata); m_free(roc);
    v_free(v);
    free(fpr); free(fnr);
    fclose(fp);
}

void face_train(gauss_t *color, char *ifname, char *rfname) {
    size_t i, s;
    double d, t;
    char fnbuf[32];
    Image *img, *ref, *out;
    MAT *tdata;
    VEC *v;
    
    // load input and refrence-mask
    img = load_image(ifname);
    ref = load_image(rfname);

    // mask input for MLE
    out = image_and(img, ref, NULL);

    // write masked image
    sprintf(fnbuf, "%s_masked.ppm", ifname);
    write_image(fnbuf, out);

    printf("Getting '%s' RG vectors for training...\n", rfname);
    // get masked RG dataset
    color->dataset = image_rgmat(out, color->dataset);

    // trim zero data
    trim_zeros(color->dataset);

    printf("Estimating color distribution for %llu samples...\n", color->dataset->m);
    // estimate distribution over masked dataset
    compute_mle(color, color->dataset->m);

    del_image(img); del_image(ref); del_image(out);
    m_free(tdata);
    v_free(v);
}

void face_exp1() {
    gauss_t color = {
        .id = 420,
        .mu = v_get(2),
        .sigma = m_get(2, 2),
        .dataset = MNULL
    };
    double c;

    face_train(&color, "train1.ppm", "ref1.ppm");

    c = 1 / (2 * M_PI * sqrt(m_det(color.sigma)));

    face_detect(&color, "train3.ppm", "ref3.ppm", 20, c / 20);
    face_detect(&color, "train6.ppm", "ref6.ppm", 20, c / 20);
}

int main(void) {
    int i, k;
    FILE *da;
    FILE *db;
    MAT *mle = MNULL;
    MAT *adata = m_get(ADATA_LEN, 2);
    MAT *bdata = m_get(BDATA_LEN, 2);
    MAT *datasets[] = { adata, bdata };
    gauss_t c1, c2;
    gauss_t *classes[] = { &c1, &c2 };
    batch_t b1 = {
        .c = 2,
        .d = 2,
        .g = case3_disc,
        .dist = classes
    };
    
    c1 = (gauss_t) {
        .id = 1,
        .mu = v_get(2),
        .sigma = m_get(2, 2),
        .dataset = adata,
        .prior = 0.3
    };

    c2 = (gauss_t) {
        .id = 2,
        .mu = v_get(2),
        .sigma = m_get(2, 2),
        .dataset = bdata,
        .prior = 0.7
    };

    c1.mu->ve[0] = 1; c1.mu->ve[1] = 1;
    m_set_val(c1.sigma, 0, 0, 1); m_set_val(c1.sigma, 0, 1, 0);
    m_set_val(c1.sigma, 1, 0, 0); m_set_val(c1.sigma, 1, 1, 1);

    c2.mu->ve[0] = 4; c2.mu->ve[1] = 4;
    m_set_val(c2.sigma, 0, 0, 1); m_set_val(c2.sigma, 0, 1, 0);
    m_set_val(c2.sigma, 1, 0, 0); m_set_val(c2.sigma, 1, 1, 1);

    da = setup_dataset("data_1A.mat", &c1, ADATA_LEN);
    db = setup_dataset("data_1B.mat", &c2, BDATA_LEN);
    if (!da || !db) {
        printf("Exiting...\n");
        return 1;
    }

    // exp 1
    // printf("============ <EXPERIMENT 1> ============\n\n");
    // mle_exp(datasets, &b1, 1);
    // printf("\n============ </EXPERIMENT 1> ============\n\n");

    c2.mu->ve[0] = 4; c2.mu->ve[1] = 4;
    m_set_val(c2.sigma, 0, 0, 4); m_set_val(c2.sigma, 0, 1, 0);
    m_set_val(c2.sigma, 1, 0, 0); m_set_val(c2.sigma, 1, 1, 8);

    da = setup_dataset("data_2A.mat", &c1, ADATA_LEN);
    db = setup_dataset("data_2B.mat", &c2, BDATA_LEN);
    if (!da || !db) {
        printf("Exiting...\n");
        return 1;
    }

    // exp 2
    // printf("============ <EXPERIMENT 2> ============\n\n");
    // mle_exp(datasets, &b1, 2);
    // printf("\n============ </EXPERIMENT 2> ============\n\n");

    // exp 3
    face_exp1();

    fclose(da); fclose(db);

    return 0;
}

void mle_exp(MAT **datasets, batch_t *batch, int it) {
    size_t t, c, s, l;
    int i, k, j;
    double pct;
    char fnbuf[32], tbuf[48];

    printf("Plotting datasets...\n");
    for (i = 0; i < batch->c; i += 1) {
        printf("\tDataset %c -> ", 'A' + i);
        sprintf(fnbuf, "plot_%d%c", it, 'A' + i);
        sprintf(tbuf, "DATASET %c - CLASS %d", 'A' + i, batch->dist[i]->id);
        plot_data(&datasets[i], 1, fnbuf, tbuf);
        printf("Done\n");
    }

    printf("\tAll Datasets -> ");
    sprintf(fnbuf, "plot_%dALL", it);
    sprintf(tbuf, "DATASETS A thru %c - ALL CLASSES", 'A' + batch->c - 1);
    plot_data(datasets, batch->c, fnbuf, tbuf);
    printf("Done\n");

    // ML estimation for exp 1 w/ differing amounts of data
    for (k = 0, pct = 1; k <= 5; k += 1) {
        if (k) {
            printf("\n============ <MLE ESTIMATES %.02lf%% data> ============\n\n", pct * 100);

            // compute ML estimations over each class and associated dataset
            for (i = 0; i < batch->c; i += 1) {
                l = pct * batch->dist[i]->dataset->m;
                printf("Dataset %d: %llu samples\n", batch->dist[i]->id, l);

                compute_mle(batch->dist[i], l);

                printf("\n~~~ MEAN (c = %d) ~~~\n", batch->dist[i]->id);
                v_output(batch->dist[i]->mu);
                printf("~~~ COVARIANCE (c = %d) ~~~\n", batch->dist[i]->id);
                m_output(batch->dist[i]->sigma);
                printf("\n");
            }

            pct *= 0.1;
        }

        // classify each class' dataset with estimated dist-parameters
        for (i = 0, t = 0, s = 0; i < batch->c; i += 1) {
            sprintf(batch->bname, "batch_%d%c-%d", it, 'A' + i, k);
            batch->n = datasets[i]->m;
            c = classify(batch, datasets[i], batch->dist[i]->id);
            t += c;
            s += datasets[i]->m;
            plot_batch(batch);
        }

        printf("\nCorrectly classified %llu of %llu (%lf%%) in total\n", t, s, 100 * (double)t / s);
    }
}

size_t classify(batch_t *batch, MAT *data, int class) {
    double *p = calloc(sizeof(double), data->n);
    size_t i, k, pm, ct;
    char fname[32];
    VEC *v;

    v = v_get(data->n);
    sprintf(fname, "%s%s.dat", DATA_DIR, batch->bname);
    batch->fdata = fopen(fname, "w+");

    for (i = 0, ct = 0; i < batch->n; i += 1) {
        // get data sample
        get_row(data, i, v);
        for (k = 0, pm = 0; k < batch->c; k += 1) {
            // compute likelihood against each class
            p[k] = batch->g(v, batch->dist[k]);
            // keep track of greatest likelihood
            if (p[k] <= p[pm]) pm = k;
        }

        // record sample classification
        for (k = 0; k < v->dim; k += 1) {
            fprintf(batch->fdata, "%lf ", v->ve[k]);
        }
        fprintf(batch->fdata, "%d\n", batch->dist[pm]->id);

        // increment correct-count if max. likelihood == correct-class id
        if (batch->dist[pm]->id == class) ct += 1; 
    }

    fclose(batch->fdata);
    v_free(v);
    free(p);

    printf("%s correctly classified %llu of %llu (%lf%%)\n", batch->bname, ct, batch->n, 100 * (double)ct / batch->n);

    return ct;
}

void plot_data(MAT **data, int n, char *fname, char *title) {
    gnuplot_ctrl *plot;
    size_t i, k, l;
    char pbuf[32], *p;

    plot = gnuplot_init();

    gnuplot_cmd(plot, "set terminal pngcairo size 800,600 enhanced font 'Consolas,12'");
    gnuplot_cmd(plot, "set output '%s%s.png'", PLOT_DIR, fname);

    gnuplot_set_xlabel(plot, "X");
    gnuplot_set_ylabel(plot, "Y");
    gnuplot_cmd(plot, "set title \"%s\"", title);

    qsort(data, n, sizeof(MAT *), dataset_len_cmp);
    
    for (l = 0; l < n; l += 1) {
        gnuplot_cmd(plot, "%s'-' using 1:%d title \"%d\" with points pointtype %d pointsize 1%s", \
                        (l == 0) ? "plot " : " ",
                        data[l]->n,
                        l + 1,
                        l + 4,
                        (n > 1 && l < n - 1) ? ",\\" : "\n");
    }


    for (l = 0; l < n; l += 1) {
        for (i = 0; i < data[l]->m; i += 1) {
            for (k = 0, p = pbuf; k < data[l]->n; k += 1) {
                p += sprintf(p, "%lf ", data[l]->me[i][k]);
            }
            gnuplot_cmd(plot, "%s", pbuf);
        }
        gnuplot_cmd(plot, "e\n");
    }

    gnuplot_close(plot);
}

void plot_batch(batch_t *b) {
    gnuplot_ctrl *plot;

    plot = gnuplot_init();
    gnuplot_cmd(plot, "set terminal pngcairo size 800,600 enhanced font 'Consolas,12'");
    gnuplot_cmd(plot, "set output '%s%s.png'", PLOT_DIR, b->bname);

    gnuplot_set_xlabel(plot, "X");
    gnuplot_set_ylabel(plot, "Y");
    gnuplot_cmd(plot, "set title \"%s\"", b->bname);


    gnuplot_cmd(plot, "plot for [i=1:%d] '%s%s.dat' using 1:%d:(column(%d) == i ? i : NaN) with \
                        points pointtype i+1 pointsize 2 lc var title sprintf(\"Class %%d\", i)", \
                        b->d, 
                        DATA_DIR, b->bname,
                        b->d, 
                        b->d + 1);

    gnuplot_close(plot);
}
