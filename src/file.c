/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <stdlib.h> /* malloc */
#include <stdio.h> /* printf */
#include <stddef.h> /* offsetof */
#include <stdint.h> /* SIZE_MAX */
#include <string.h> /* strlen, strncat */
#include <ctype.h> /* isspace */
#include <zlib.h> /* gzopen, gzgets, gzeof, gzclose */

#include "file.h"

#ifdef UNIT_TESTING
extern void * _test_malloc(const size_t size, const char * file, const int line);
#define malloc(size) _test_malloc(size, __FILE__, __LINE__)
#endif

#define MAXCHAR 1024
void gzreadnext(gzFile F, char data[], int n)
{
    char * ptr = gzgets(F, data, n);
    if(ptr == Z_NULL) {
        data[0] = '\0';
        return;
    }
    int cnt;
    for(cnt = 0; (data[cnt] != '\0') && (data[cnt] != '\n') && (cnt < n - 1); cnt++) {
        /* netgen is inconsistent in its use of \r and \n */
        if (data[cnt] == '\r') {
            data[cnt] = '\n';
            break;
        }
    }
    data[cnt + 1] = '\0';
}

#define LEN 32
enum fileformat_req {OPTIONAL, REQUIRED, REQUIRED_FIRST};
typedef struct {
    char section [LEN];
    int n_size;
    size_t offset;
    int (*sscanf_funcptr)(const char *, model *, mesh *, int i);
    int (*test_funcptr)(model *, mesh *);
    enum fileformat_req req;
    int cnt;
} fileformat;

int check_req(const fileformat * f, enum fileformat_req needs)
{
    int i;
    for(i = 0; f[i].section[0]; i++) {
        if((f[i].req >= needs) && (!f[i].cnt)) {
            return i + 1;
        }
    }
    return 0;
}

int readfile_loop(char filename[], model * m, fileformat * f_list)
{
    if(m == NULL) {
        return 1;
    }
    int ret = 1;
    int cnt;
    char data[MAXCHAR];
    model_free(m);
    gzFile F = gzopen(filename, "r");
    if(!F) {
        printf("err: failed to open %s\n", filename);
        goto __quit;
    }
    printf("reading %s\n", filename);
    mesh * mm = &(m->fwd);
    int first = 1;
    do {
        gzreadnext(F, data, MAXCHAR);
        fileformat * f = &(f_list[0]);
        while(f->section[0]) {
            char section_n [LEN + 1];
            char * eol = stpncpy(section_n, f->section, LEN + 1);
            *eol = '\n';
            if(strcmp(data, section_n) != 0) {
                if(first && (f->req == REQUIRED_FIRST)) {
                    printf("err: %s: missing header %s\n", filename, f->section);
                    goto __quit;
                }
                f++;
                continue;
            }
            first = 0;
            f->cnt += 1;
            if(!f->n_size) {
                printf("%s\n", f->section);
                f++;
                break;
            }
            gzreadnext(F, data, MAXCHAR);
            int * n = (int *) ((char *)mm + f->offset);
            cnt = 0;
            int idx = 0;
            int j;
            for(j = 0; j < f->n_size; j++) {
                cnt += sscanf(&(data[idx]), "%d %n", &(n[j]), &idx);
            }
            if(cnt != f->n_size) {
                printf("err: bad %s\n", f->section);
                goto __quit;
            }
            printf("%s %d\n", f->section, *n);
            if(f->sscanf_funcptr) {
                int i;
                int * n = (int *) ((char *)mm + f->offset);
                for(i = 0; i < *n; i += 1) {
                    gzreadnext(F, data, MAXCHAR);
                    int ret = (*f->sscanf_funcptr)(data, m, mm, i);
                    if( ret ) {
                        if(ret == 1) {
                            printf("err: bad %s\n", f->section);
                        }
                        else {
                            printf("err: out of memory %s\n", f->section);
                        }
                        goto __quit;
                    }
                }
            }
            if(f->test_funcptr && (*f->test_funcptr)(m, mm)) {
                printf("err: bad %s\n", f->section);
                goto __quit;
            }
        }
    } while(!gzeof(F));
    int chk = check_req(f_list, REQUIRED);
    if(chk) {
        printf("err: %s: missing %s\n", filename, (f_list[chk - 1]).section);
        goto __quit;
    }
    if( (mm->n_nodes == 0) || (mm->n_elems == 0) || (mm->n_se == 0) ) {
        printf("err: %s: empty mesh nodes = %d, elems = %d, surfs = %d\n",
               filename, mm->n_nodes, mm->n_elems, mm->n_se);
        ret = 2;
    }
    else {
        ret = 0;
    }
__quit:
    if(ret != 0) {
        model_free(m);
    }
    if(F) {
        if(gzclose(F) != Z_OK) {
            printf("err: failed to close %s\n", filename);
            ret = 1;
        }
    }
    return ret;
}

int readngvol_test_dimension(model * m, mesh * mm)
{
    return !((mm->dim == 3) || (mm->dim == 2));
}

int readngvol_test_geomtype(model * m, mesh * mm)
{
    return (mm->type != 0);
}

int readngvol_surfaceelements(const char * data, model * m, mesh * mm, const int se)
{
    if(!mm->surfaceelems) {
        mm->surfaceelems = malloc(sizeof(int) * mm->n_se * mm->dim);
    }
    if(!mm->bc) {
        mm->bc = malloc(sizeof(int) * mm->n_se * 1);
    }
    if(!mm->surfaceelems || !mm->bc) {
        return 2;
    }
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &mm->bc[se], &mm->surfaceelems[3 * se + 0], &mm->surfaceelems[3 * se + 1], &mm->surfaceelems[3 * se + 2]);
    }
    else {
        cnt = sscanf(data, "%*d %d %*d %*d 2 %d %d\n", &mm->bc[se], &mm->surfaceelems[2 * se + 0], &mm->surfaceelems[2 * se + 1]);
    }
    return (cnt != mm->dim + 1);
}

int readngvol_points(const char * data, model * m, mesh * mm, const int node)
{
    if(!mm->nodes) {
        mm->nodes = malloc(sizeof(double) * mm->n_nodes * mm->dim);
    }
    if(!mm->nodes) {
        return 2;
    }
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, " %lf %lf %lf\n", &mm->nodes[3 * node + 0], &mm->nodes[3 * node + 1], &mm->nodes[3 * node + 2]);
    }
    else {
        cnt = sscanf(data, " %lf %lf\n", &mm->nodes[2 * node + 0], &mm->nodes[2 * node + 1]);
    }
    return (cnt != mm->dim);
}

int readngvol_volumeelements(const char * data, model * m, mesh * mm, const int elem)
{
    if(!mm->elems) {
        mm->elems = malloc(sizeof(int) * mm->n_elems * (mm->dim + 1));
    }
    if(!mm->matidx) {
        mm->matidx = malloc(sizeof(int) * mm->n_elems * 1);
    }
    if(!mm->elems || !mm->matidx) {
        return 2;
    }
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, "%d 4 %d %d %d %d\n", &mm->matidx[elem], &mm->elems[4 * elem + 0], &mm->elems[4 * elem + 1], &mm->elems[4 * elem + 2], &mm->elems[4 * elem + 3]);
    }
    else {
        cnt = sscanf(data, "%d 3 %d %d %d\n", &mm->matidx[elem], &mm->elems[3 * elem + 0], &mm->elems[3 * elem + 1], &mm->elems[3 * elem + 2]);
    }
    return (cnt != mm->dim + 2);
}

int readzh_test_format1(model * m, mesh * mm)
{
    return (m->format != 1);
}

int readzh_stimmeas(const char * data, model * m, mesh * mm, const int i)
{
    if(!m->stimmeas) {
        m->stimmeas = malloc(sizeof(int) * m->n_stimmeas * 4);
    }
    if(!m->stimmeas) {
        return 2;
    }
    int cnt = sscanf(data, "%d %d %d %d\n",
                     &m->stimmeas[4 * i + 0], &m->stimmeas[4 * i + 1],
                     &m->stimmeas[4 * i + 2], &m->stimmeas[4 * i + 3]);
    return (cnt != 4);
}

int readzh_data(const char * data, model * m, mesh * mm, const int i)
{
    const int rows = m->n_data[0];
    const int cols = m->n_data[1];
    if(!m->data) {
        m->data = malloc(sizeof(double) * rows * cols);
    }
    if(!m->data) {
        return 2;
    }
    int cnt = 0;
    int idx = 0;
    int j;
    for(j = 0; j < cols; j++) {
        cnt += sscanf(&(data[idx]), " %lf %n", &(m->data[cols * i + j]), &idx);
    }
    return (cnt != cols);
}

int readzh_params(const char * data, model * m, mesh * mm, const int i)
{
    const int rows = m->n_params[0];
    const int cols = m->n_params[1];
    if(!m->params) {
        m->params = malloc(sizeof(double) * rows * cols);
    }
    if(!m->params) {
        return 2;
    }
    int cnt = 0;
    int idx = 0;
    int j;
    for(j = 0; j < cols; j++) {
        cnt += sscanf(&(data[idx]), " %lf %n", &(m->params[cols * i + j]), &idx);
    }
    return (cnt != cols);
}

int readfile_ngvol(char filename[], model * m)
{
    fileformat ngvol [] = {
        {"mesh3d", 0, 0, NULL, NULL, REQUIRED_FIRST},
        {"dimension", 1, offsetof(model, fwd.dim), NULL, &readngvol_test_dimension, REQUIRED},
        {"geomtype", 1, offsetof(model, fwd.type), NULL, &readngvol_test_geomtype, REQUIRED},
        /* GEOM_STL=11: surfaceelementsgi
         * GEOM_OCC=12, GEOM_ACIS=13: surfaceelementsuv
         * default: surfaceelements
         */
        {"surfaceelements", 1, offsetof(model, fwd.n_se), &readngvol_surfaceelements, NULL, REQUIRED},
        {"points", 1, offsetof(model, fwd.n_nodes), &readngvol_points, NULL, REQUIRED},
        {"volumeelements", 1, offsetof(model, fwd.n_elems), &readngvol_volumeelements, NULL, REQUIRED},
        {{0}}
    };
    return readfile_loop(filename, m, ngvol);
}

int readfile(char filename[], model * m)
{
    fileformat zh_format1 [] = {
        {"zedhat", 0, 0, NULL, NULL, REQUIRED_FIRST},
        {"format", 1, offsetof(model, format), NULL, &readzh_test_format1, REQUIRED},
        /* netgen .vol format */
        {"dimension", 1, offsetof(model, fwd.dim), NULL, &readngvol_test_dimension, REQUIRED},
        {"geomtype", 1, offsetof(model, fwd.type), NULL, &readngvol_test_geomtype, REQUIRED},
        {"surfaceelements", 1, offsetof(model, fwd.n_se), &readngvol_surfaceelements, NULL, REQUIRED},
        {"points", 1, offsetof(model, fwd.n_nodes), &readngvol_points, NULL, REQUIRED},
        {"volumeelements", 1, offsetof(model, fwd.n_elems), &readngvol_volumeelements, NULL, REQUIRED},
        {"stimmeas", 1, offsetof(model, n_stimmeas), &readzh_stimmeas, NULL, OPTIONAL},
        {"data", 2, offsetof(model, n_data), &readzh_data, NULL, OPTIONAL},
        {"parameters", 2, offsetof(model, n_params), &readzh_params, NULL, OPTIONAL},
        {{0}}
    };
    return readfile_loop(filename, m, zh_format1);
}
