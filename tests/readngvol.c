
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf, feof, fgets */
#include <libgen.h> /* basename */
#include <string.h> /* strlen */
#include <ctype.h> /* isspace */

#define MAXCHAR 1024
int main(int argc, char ** argv)
{
    int ret = 0;
    int cnt;
    char data[MAXCHAR];
    double * nodes = NULL; int n_nodes;
    int * elems = NULL; int * matidx = NULL; int n_elems;
    int * surfaceelems = NULL; int * bc = NULL; int n_se;
    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol>\n", basename(argv[0]));
        return 1;
    }
    FILE * F = fopen(argv[1], "r");
    if(!F) {
        fprintf(stderr, "error: failed to open %s\n", argv[0]);
        goto __quit;
    }
    printf("reading %s\n", argv[1]);
    if(!fgets(data, 8, F)) {
        goto __quit;
    }
    /* DANGER: inconsistent use of \r and \n by netgen */
    for(cnt = 0; data[cnt] != 0; cnt++) {
        if (data[cnt] == '\r' || data[cnt] == '\n') {
            data[cnt] = '\0';
        }
    }
    if(strcmp(data, "mesh3d")) {
        fprintf(stderr, "err: bad header\n");
        goto __quit;
    }
    while(!feof(F)) {
        if(!fgets(data, MAXCHAR, F)) {
            break;
        }
        for(cnt = 0; data[cnt] != '\0'; cnt++) {
            if (data[cnt] == '\r' || data[cnt] == '\n') {
                data[cnt] = '\0';
            }
        }
        if(!strcmp(data, "dimension")) {
            int dim = 0;
            cnt = fscanf(F, "%d\n", &dim);
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
        else if(!strcmp(data, "geomtype")) {
            int type = -1;
            cnt = fscanf(F, "%d\n", &type);
            if((cnt != 1) || (type != 0)) {
                fprintf(stderr, "err: bad geomtype\n");
                goto __quit;
            }
            printf("geomtype %d\n", type);
        }
        else if(!strcmp(data, "surfaceelements")) {
            cnt = fscanf(F, "%d\n", &n_se);
            printf("%d surfaceelements\n", n_se);
            int se = 0;
            surfaceelems = malloc(sizeof(int) * n_se * 3);
            bc = malloc(sizeof(int) * n_se * 1);
            if(!surfaceelems || !bc) {
                return -1;
            }
            while(!feof(F) && se < n_se ) {
                cnt += fscanf(F, "%*d %d %*d %*d 3 %d %d %d\n", &bc[se], &surfaceelems[3 * se + 0], &surfaceelems[3 * se + 1], &surfaceelems[3 * se + 2]);
                se++;
            }
            if(cnt != n_se * 4 + 1) {
                fprintf(stderr, "err: bad surfaceelements\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "points")) {
            cnt = fscanf(F, "%d\n", &n_nodes);
            printf("%d points\n", n_nodes);
            int node = 0;
            nodes = malloc(sizeof(double) * n_nodes * 3);
            if(!nodes) {
                return -1;
            }
            while(!feof(F) && node < n_nodes ) {
                cnt += fscanf(F, " %lf %lf %lf\n", &nodes[3 * node + 0], &nodes[3 * node + 1], &nodes[3 * node + 2]);
                node++;
            }
            if(cnt != n_nodes * 3 + 1) {
                fprintf(stderr, "err: bad points\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "volumeelements")) {
            cnt = fscanf(F, "%d\n", &n_elems);
            printf("%d volumeelements\n", n_elems);
            int elem = 0;
            elems = malloc(sizeof(int) * n_elems * 4);
            matidx = malloc(sizeof(int) * n_elems * 1);
            if(!elems || !matidx) {
                return -1;
            }
            while(!feof(F) && elem < n_elems ) {
                cnt += fscanf(F, "%d 4 %d %d %d %d\n", &matidx[elem], &elems[4 * elem + 0], &elems[4 * elem + 1], &elems[4 * elem + 2], &elems[4 * elem + 3]);
                elem++;
            }
            if(cnt != n_elems * 5 + 1) {
                fprintf(stderr, "err: bad elems\n");
                goto __quit;
            }
        }
    }
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
        ret = fclose(F);
        if(ret) {
            fprintf(stderr, "error: failed to close %s\n", argv[1]);
        }
    }
    return ret;
}
