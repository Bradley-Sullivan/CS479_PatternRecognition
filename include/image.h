#ifndef IMAGE_H
#define IMAGE_H

#define IMAGE_DIR "images/"
#define MAX_FPATH 256

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>

typedef struct Image {
    uint16_t m, n, q;
    size_t size;
    uint8_t *data;
} Image;

Image *image_and(Image *a, Image *b, Image *out);
Image *new_image(int m, int n, int q);
Image *copy_image(Image* img);
int del_image(Image *img);
Image *load_image(const char *fname);
int load_header(FILE *fp, Image *img);
int load_data(FILE *fp, Image *img);
int write_image(const char *fname, Image *img);

#endif // IMAGE_H