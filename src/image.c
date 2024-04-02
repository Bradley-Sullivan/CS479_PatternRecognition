#include "image.h"

Image *image_and(Image *a, Image *b, Image *out) {
    Image *result;

    if (!out) result = new_image(a->m, a->n, a->q);
    else result = out;

    for (size_t i = 0; i < a->size; i += 1) {
        result->data[i] = (uint8_t) a->data[i] & b->data[i];
    }
    
    return result;
}

Image *new_image(int m, int n, int q) {
    Image *ret = (Image*)malloc(sizeof(Image));

    ret->m = m;
    ret->n = n;
    ret->q = q;
    ret->size = m * n * 3;

    ret->data = (uint8_t*) calloc(ret->size, sizeof(uint8_t));

    return ret;
}

Image *copy_image(Image* img) {
    Image *ret = new_image(img->m, img->n, img->q);

    if (!memcpy(ret->data, img->data, img->size)) {
        fprintf(stderr, "Error. Failed to copy input data array.\n");
        return NULL;
    }

    return ret;
}

int del_image(Image *img) {
    if (!img) {
        fprintf(stderr, "Cannot free NULL Image pointer.\n");
        return 1;
    } else if (!img->data) {
        fprintf(stderr, "Cannot free NULL Image-data pointer.\n");
        return 1;
    }

    free(img->data);

    free(img);

    return 0;
}

Image *load_image(const char *fname) {
    char fpath[64];    

    sprintf(fpath, "%s%s", IMAGE_DIR, fname);

    Image *img = (Image*) malloc(sizeof(Image));
    FILE *fp = fopen(fpath, "r");

    if (!fp) { 
        fprintf(stderr, "Error opening '%s'.\n", fname); 
        return NULL;
    }

    // gets file size
    fseek(fp, 0, SEEK_END);
    size_t fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    uint8_t *fbuf = (uint8_t*) malloc(sizeof(uint8_t) * fsize);

    // loads entire file into memory
    if (!fread(fbuf, sizeof(uint8_t), fsize, fp)) {
        fprintf(stderr, "Error loading file into buffer.\n");
        return NULL;
    }

    // opens filestream on file buffer
    FILE *fbuf_stream = fmemopen(fbuf, fsize, "r");

    // extracts header information from file buffer
    if (load_header(fbuf_stream, img)) {
        fprintf(stderr, "Error loading file header.\n");
        return NULL;
    }

    if (load_data(fbuf_stream, img)) {
        fprintf(stderr, "Error loading file data.\n");
        return NULL;
    }
    
    free(fbuf);
    fclose(fp);

    return img;
}

int load_header(FILE *fp, Image *img) {
    char *rbuf = (char*) malloc(sizeof(char) * 64);
    int rgb = 0;

    fgets(rbuf, 64, fp);

    if (strcmp(rbuf, "P6\n") == 0) {
        rgb = 1;
    }

    do { fgets(rbuf, 64, fp); } while (rbuf[0] == '#');

    sscanf(rbuf, "%hu %hu\n", &img->m, &img->n);

    fscanf(fp, "%hu\n", &img->q);

    img->size = img->m * img->n * ((rgb) ? 3 : 1);

    free(rbuf);

    return 0;
}

int load_data(FILE *fp, Image *img) {
    uint8_t *img_data = (uint8_t*) malloc(sizeof(uint8_t) * img->size);

    if (!fread(img_data, sizeof(uint8_t), img->size, fp)) {
        fprintf(stderr, "Error. Image is not %d by %d.\n", img->m, img->n);
        return 1;
    }

    img->data = img_data;

    return 0;
}

int write_image(const char *fname, Image *img) {
    char pbuf[MAX_FPATH];
    FILE *fp;

    sprintf(pbuf, "%s%s", IMAGE_DIR, fname);
    
    fp = fopen(pbuf, "w");
    if (!fp) { fprintf(stderr, "Error opening destination file '%s'.\n", fname); return 1; }

    // writes header
    fprintf(fp, "P%d\n%d %d\n%d\n", (img->size > img->m * img->n) ? 6 : 5, img->m, img->n, img->q);

    // writes data
    size_t wsize = fwrite(img->data, sizeof(uint8_t), img->size, fp);
    if (wsize != img->size) {
        fprintf(stderr, "Error writing image data to file.\n");
        return 1;
    }

    fclose(fp);

    return 0;
}