
#include <stdlib.h> /* malloc */
#include <stdio.h> /* printf */
#include <string.h> /* strlen */
#include <ctype.h> /* isspace */
#include <zlib.h> /* gzopen, gzgets, gzeof, gzclose */

#include "readngvol.h"

#ifdef UNIT_TESTING
extern void * _test_malloc(const size_t size, const char * file, const int line);
#define malloc(size) _test_malloc(size, __FILE__, __LINE__)
#endif


#define MAXCHAR 1024
void gzreadnext(gzFile F, char data[], int n)
{
    data[0] = '\0';
    gzgets(F, data, n);
    /* DANGER: inconsistent use of \r and \n by netgen */
    int cnt;
    for(cnt = 0; (data[cnt] != '\0') && (data[cnt] != '\n') && (cnt < n); cnt++) {
        if (data[cnt] == '\r') {
            data[cnt] = '\n';
            break;
        }
    }
    data[cnt + 1] = '\0';
}

int readngvol(char filename[], mesh * m)
{
    int ret = 1;
    int cnt;
    char data[MAXCHAR];
    mesh_free(m);
    gzFile F = gzopen(filename, "r");
    if(!F) {
        printf("err: failed to open %s\n", filename);
        goto __quit;
    }
    printf("reading %s\n", filename);
    gzreadnext(F, data, 8);
    if(strcmp(data, "mesh3d\n")) {
        printf("err: bad header\n");
        goto __quit;
    }
    while(!gzeof(F)) {
        gzreadnext(F, data, MAXCHAR);
        if(!strcmp(data, "dimension\n")) {
            gzreadnext(F, data, MAXCHAR);
            cnt = sscanf(data, "%d\n", &(m->dim));
            if((cnt != 1) || (m->dim != 3)) {
                printf("err: bad dimension\n");
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
            gzreadnext(F, data, MAXCHAR);
            cnt = sscanf(data, "%d\n", &type);
            if((cnt != 1) || (type != 0)) {
                printf("err: bad geomtype\n");
                goto __quit;
            }
            printf("geomtype %d\n", type);
        }
        else if(!strcmp(data, "surfaceelements\n")) {
            gzreadnext(F, data, MAXCHAR);
            cnt = sscanf(data, "%d\n", &m->n_se);
            printf("%d surfaceelements\n", m->n_se);
            int se = 0;
            m->surfaceelems = malloc(sizeof(int) * m->n_se * 3);
            m->bc = malloc(sizeof(int) * m->n_se * 1);
            if(!m->surfaceelems || !m->bc) {
                goto __quit;
            }
            while(!gzeof(F) && se < m->n_se ) {
                gzreadnext(F, data, MAXCHAR);
                cnt += sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &m->bc[se], &m->surfaceelems[3 * se + 0], &m->surfaceelems[3 * se + 1], &m->surfaceelems[3 * se + 2]);
                se++;
            }
            if(cnt != m->n_se * 4 + 1) {
                printf("err: bad surfaceelements\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "points\n")) {
            gzreadnext(F, data, MAXCHAR);
            cnt = sscanf(data, "%d\n", &m->n_nodes);
            printf("%d points\n", m->n_nodes);
            int node = 0;
            m->nodes = malloc(sizeof(double) * m->n_nodes * 3);
            if(!m->nodes) {
                goto __quit;
            }
            while(!gzeof(F) && node < m->n_nodes ) {
                gzreadnext(F, data, MAXCHAR);
                cnt += sscanf(data, " %lf %lf %lf\n", &m->nodes[3 * node + 0], &m->nodes[3 * node + 1], &m->nodes[3 * node + 2]);
                node++;
            }
            if(cnt != m->n_nodes * 3 + 1) {
                printf("err: bad points\n");
                goto __quit;
            }
        }
        else if(!strcmp(data, "volumeelements\n")) {
            gzreadnext(F, data, MAXCHAR);
            cnt = sscanf(data, "%d\n", &m->n_elems);
            printf("%d volumeelements\n", m->n_elems);
            int elem = 0;
            m->elems = malloc(sizeof(int) * m->n_elems * 4);
            m->matidx = malloc(sizeof(int) * m->n_elems * 1);
            if(!m->elems || !m->matidx) {
                goto __quit;
            }
            while(!gzeof(F) && elem < m->n_elems ) {
                gzreadnext(F, data, MAXCHAR);
                cnt += sscanf(data, "%d 4 %d %d %d %d\n", &m->matidx[elem], &m->elems[4 * elem + 0], &m->elems[4 * elem + 1], &m->elems[4 * elem + 2], &m->elems[4 * elem + 3]);
                elem++;
            }
            if(cnt != m->n_elems * 5 + 1) {
                printf("err: bad elems\n");
                goto __quit;
            }
        }
    }
    if( (m->n_nodes == 0) || (m->n_elems == 0) || (m->n_se == 0) ) {
        printf("err: %s: empty mesh nodes = %d, elems = %d, surfs = %d\n",
                filename, m->n_nodes, m->n_elems, m->n_se);
        ret = 2;
    }
    else {
        ret = 0;
    }
__quit:
    if(ret != 0) {
        mesh_free(m);
    }
    if(F) {
        if(gzclose(F) != Z_OK) {
            printf("err: failed to close %s\n", filename);
            ret = 1;
        }
    }
    return ret;
}
