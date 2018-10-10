
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf */
#include <libgen.h> /* basename */
#include <string.h> /* strlen */
#include <ctype.h> /* isspace */
#include <zlib.h> /* gzopen, gzgets, gzeof, gzclose */

#define MAXCHAR 1024
int gzreadnext(gzFile F, char data[], int n)
{
    if(!gzgets(F, data, n)) {
        return 0; /* error */
    }
    /* DANGER: inconsistent use of \r and \n by netgen */
    int cnt;
    for(cnt = 0; (data[cnt] != '\0') && (data[cnt] != '\n') && (cnt < n); cnt++) {
        if (data[cnt] == '\r') {
            data[cnt] = '\n';
            break;
        }
    }
    data[cnt + 1] = '\0';
    return 1; /* success */
}

int main(int argc, char ** argv)
{
    int ret = 1;
    int cnt;
    char data[MAXCHAR];
    double * nodes = NULL; int n_nodes;
    int * elems = NULL; int * matidx = NULL; int n_elems;
    int * surfaceelems = NULL; int * bc = NULL; int n_se;
    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol>\n", basename(argv[0]));
        return 1;
    }
    gzFile F = gzopen(argv[1], "r");
    if(!F) {
        fprintf(stderr, "error: failed to open %s\n", argv[1]);
        goto __quit;
    }
    printf("reading %s\n", argv[1]);
    if(!gzreadnext(F, data, 8) || strcmp(data, "mesh3d\n")) {
        fprintf(stderr, "err: bad header\n");
        goto __quit;
    }
    while(!gzeof(F)) {
        if(!gzreadnext(F, data, MAXCHAR)) {
            goto __quit;
        }
        if(!strcmp(data, "dimension\n")) {
            int dim = 0;
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &dim);
            if((cnt != 1) || (dim != 3)) {
                fprintf(stderr, "err: bad dimension\n");
                goto __quit;
            }
            printf("dimension %d\n", dim);
        }
        /* GEOM_STL=11: surfaceelementsgi
         * GEOM_OCC=12, GEOM_ACIS=13: surfaceelementsuv
         * default: surfaceelements
         */
        else if(!strcmp(data, "geomtype\n")) {
            int type = -1;
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &type);
            if((cnt != 1) || (type != 0)) {
                fprintf(stderr, "err: bad geomtype\n");
                goto __quit;
            }
            printf("geomtype %d\n", type);
        }
        else if(!strcmp(data, "surfaceelements\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &n_se);
            printf("%d surfaceelements\n", n_se);
            int se = 0;
            surfaceelems = malloc(sizeof(int) * n_se * 3);
            bc = malloc(sizeof(int) * n_se * 1);
            if(!surfaceelems || !bc) {
                goto __quit;
            }
            while(!gzeof(F) && se < n_se ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &bc[se], &surfaceelems[3 * se + 0], &surfaceelems[3 * se + 1], &surfaceelems[3 * se + 2]);
                se++;
            }
            if(cnt != n_se * 4 + 1) {
                fprintf(stderr, "err: bad surfaceelements\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "points\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &n_nodes);
            printf("%d points\n", n_nodes);
            int node = 0;
            nodes = malloc(sizeof(double) * n_nodes * 3);
            if(!nodes) {
                return -1;
            }
            while(!gzeof(F) && node < n_nodes ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, " %lf %lf %lf\n", &nodes[3 * node + 0], &nodes[3 * node + 1], &nodes[3 * node + 2]);
                node++;
            }
            if(cnt != n_nodes * 3 + 1) {
                fprintf(stderr, "err: bad points\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "volumeelements\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &n_elems);
            printf("%d volumeelements\n", n_elems);
            int elem = 0;
            elems = malloc(sizeof(int) * n_elems * 4);
            matidx = malloc(sizeof(int) * n_elems * 1);
            if(!elems || !matidx) {
                goto __quit;
            }
            while(!gzeof(F) && elem < n_elems ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, "%d 4 %d %d %d %d\n", &matidx[elem], &elems[4 * elem + 0], &elems[4 * elem + 1], &elems[4 * elem + 2], &elems[4 * elem + 3]);
                elem++;
            }
            if(cnt != n_elems * 5 + 1) {
                fprintf(stderr, "err: bad elems\n");
                goto __quit;
            }
        }
    }
    ret = 0;
//   if((n_elems > 0) && (n_nodes > 0)) {
//      printf("happiness!\n");
//   }
__quit:
    free(nodes); nodes = NULL;
    free(elems); elems = NULL;
    free(matidx); matidx = NULL;
    free(surfaceelems); surfaceelems = NULL;
    free(bc); bc = NULL;
    if(F) {
        ret = gzclose(F);
        if(ret) {
            fprintf(stderr, "error: failed to close %s\n", argv[1]);
        }
    }
    return ret;
}
