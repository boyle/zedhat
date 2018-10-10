
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf */
#include <string.h> /* strlen */
#include <ctype.h> /* isspace */
#include <zlib.h> /* gzopen, gzgets, gzeof, gzclose */

#include "readngvol.h"

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

int readngvol(char filename[], struct mesh *m)
{
    int ret = 1;
    int cnt;
    char data[MAXCHAR];
    /* free(mesh) */
    /* init mesh */

    gzFile F = gzopen(filename, "r");
    if(!F) {
        fprintf(stderr, "error: failed to open %s\n", filename);
        goto __quit;
    }
    printf("reading %s\n", filename);

    if(!gzreadnext(F, data, 8) || strcmp(data, "mesh3d\n")) {
        fprintf(stderr, "err: bad header\n");
        goto __quit;
    }
    while(!gzeof(F)) {
        if(!gzreadnext(F, data, MAXCHAR)) {
            goto __quit;
        }
        if(!strcmp(data, "dimension\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &(m->dim));
            if((cnt != 1) || (m->dim != 3)) {
                fprintf(stderr, "err: bad dimension\n");
                goto __quit;
            }
            printf("dimension %d\n", m->dim);
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
            cnt = sscanf(data, "%d\n", &m->n_se);
            printf("%d surfaceelements\n", m->n_se);
            int se = 0;
            m->surfaceelems = malloc(sizeof(int) * m->n_se * 3);
            m->bc = malloc(sizeof(int) * m->n_se * 1);
            if(!m->surfaceelems || !m->bc) {
                goto __quit;
            }
            while(!gzeof(F) && se < m->n_se ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &m->bc[se], &m->surfaceelems[3 * se + 0], &m->surfaceelems[3 * se + 1], &m->surfaceelems[3 * se + 2]);
                se++;
            }
            if(cnt != m->n_se * 4 + 1) {
                fprintf(stderr, "err: bad surfaceelements\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "points\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &m->n_nodes);
            printf("%d points\n", m->n_nodes);
            int node = 0;
            m->nodes = malloc(sizeof(double) * m->n_nodes * 3);
            if(!m->nodes) {
                return -1;
            }
            while(!gzeof(F) && node < m->n_nodes ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, " %lf %lf %lf\n", &m->nodes[3 * node + 0], &m->nodes[3 * node + 1], &m->nodes[3 * node + 2]);
                node++;
            }
            if(cnt != m->n_nodes * 3 + 1) {
                fprintf(stderr, "err: bad points\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "volumeelements\n")) {
            if(!gzreadnext(F, data, MAXCHAR)) {
                goto __quit;
            }
            cnt = sscanf(data, "%d\n", &m->n_elems);
            printf("%d volumeelements\n", m->n_elems);
            int elem = 0;
            m->elems = malloc(sizeof(int) * m->n_elems * 4);
            m->matidx = malloc(sizeof(int) * m->n_elems * 1);
            if(!m->elems || !m->matidx) {
                goto __quit;
            }
            while(!gzeof(F) && elem < m->n_elems ) {
                if(!gzreadnext(F, data, MAXCHAR)) {
                    break;
                }
                cnt += sscanf(data, "%d 4 %d %d %d %d\n", &m->matidx[elem], &m->elems[4 * elem + 0], &m->elems[4 * elem + 1], &m->elems[4 * elem + 2], &m->elems[4 * elem + 3]);
                elem++;
            }
            if(cnt != m->n_elems * 5 + 1) {
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
    free(m->nodes); m->nodes = NULL;
    free(m->elems); m->elems = NULL;
    free(m->matidx); m->matidx = NULL;
    free(m->surfaceelems); m->surfaceelems = NULL;
    free(m->bc); m->bc = NULL;
    if(F) {
        ret = gzclose(F);
        if(ret) {
            fprintf(stderr, "error: failed to close %s\n", filename);
        }
    }
    return ret;
}
