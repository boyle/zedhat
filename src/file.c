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
    /* DANGER: inconsistent use of \r and \n by netgen */
    int cnt;
    for(cnt = 0; (data[cnt] != '\0') && (data[cnt] != '\n') && (cnt < n - 1); cnt++) {
        if (data[cnt] == '\r') {
            data[cnt] = '\n';
            break;
        }
    }
    data[cnt + 1] = '\0';
}

#define LEN 256
enum fileformat_req {OPTIONAL, REQUIRED, REQUIRED_FIRST};
typedef struct {
    char section [LEN];
    size_t offset;
    int (*sscanf)(const char *, void *, int i); /* int cnt = sscanf(const char * str, void* d) */
    int (*test)(void *); /* int fail = test(void *d) */
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

int readngvol_surfaceelements(const char * data, void * ptr, const int se);
int readngvol_volumeelements(const char * data, void * ptr, const int elem);
int readngvol_points(const char * data, void * ptr, const int node);

int readngvol_test_dimension(void * ptr)
{
    const mesh * m = ptr;
    return (m->dim != 3);
}
int readngvol_test_geomtype(void * ptr)
{
    const mesh * m = ptr;
    return (m->type != 0);
}

int readfile(char filename[], model * m)
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
    fileformat f_list [] = {
        {"mesh3d", SIZE_MAX, NULL, NULL, REQUIRED_FIRST},
        {"dimension", offsetof(mesh, dim), NULL, &readngvol_test_dimension, REQUIRED},
        {"geomtype", offsetof(mesh, type), NULL, &readngvol_test_geomtype, REQUIRED},
        /* GEOM_STL=11: surfaceelementsgi
         * GEOM_OCC=12, GEOM_ACIS=13: surfaceelementsuv
         * default: surfaceelements
         */
        {"surfaceelements", offsetof(mesh, n_se), &readngvol_surfaceelements, NULL, REQUIRED},
        {"points", offsetof(mesh, n_nodes), &readngvol_points, NULL, REQUIRED},
        {"volumeelements", offsetof(mesh, n_elems), &readngvol_volumeelements, NULL, REQUIRED},
        {{0}}
    };
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
            if(f->offset == SIZE_MAX) {
                printf("%s\n", f->section);
                f++;
                break;
            }
            gzreadnext(F, data, MAXCHAR);
            int * n = (int *) ((char *)mm + f->offset);
            cnt = sscanf(data, "%d\n", n);
            if(cnt != 1) {
                printf("err: bad %s\n", f->section);
                goto __quit;
            }
            printf("%s %d\n", f->section, *n);
            if(f->sscanf) {
                int i;
                for(i = 0; i < *n; i += 1) {
                    gzreadnext(F, data, MAXCHAR);
                    if( (*f->sscanf)(data, mm, i) ) { /* sscanf(data, mm); err=1 */
                        goto __quit;
                    }
                }
            }
            if(f->test && (*f->test)(mm)) { /* test(mm); err=1 */
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

int readngvol_surfaceelements(const char * data, void * ptr, const int se)
{
    mesh * mm = ptr;
    if(!mm->surfaceelems) {
        mm->surfaceelems = malloc(sizeof(int) * mm->n_se * 3);
    }
    if(!mm->bc) {
        mm->bc = malloc(sizeof(int) * mm->n_se * 1);
    }
    if(!mm->surfaceelems || !mm->bc) {
        return 2;
    }
    int cnt = sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &mm->bc[se], &mm->surfaceelems[3 * se + 0], &mm->surfaceelems[3 * se + 1], &mm->surfaceelems[3 * se + 2]);
    if(cnt != 4) {
        printf("err: bad surfaceelements\n");
        return 1;
    }
    return 0;
}

int readngvol_points(const char * data, void * ptr, const int node)
{
    mesh * mm = ptr;
    if(!mm->nodes) {
        mm->nodes = malloc(sizeof(double) * mm->n_nodes * 3);
    }
    if(!mm->nodes) {
        return 2;
    }
    int cnt = sscanf(data, " %lf %lf %lf\n", &mm->nodes[3 * node + 0], &mm->nodes[3 * node + 1], &mm->nodes[3 * node + 2]);
    if(cnt != 3) {
        printf("err: bad points\n");
        return 1;
    }
    return 0;
}

int readngvol_volumeelements(const char * data, void * ptr, const int elem)
{
    mesh * mm = ptr;
    if(!mm->elems) {
        mm->elems = malloc(sizeof(int) * mm->n_elems * 4);
    }
    if(!mm->matidx) {
        mm->matidx = malloc(sizeof(int) * mm->n_elems * 1);
    }
    if(!mm->elems || !mm->matidx) {
        return 2;
    }
    int cnt = sscanf(data, "%d 4 %d %d %d %d\n", &mm->matidx[elem], &mm->elems[4 * elem + 0], &mm->elems[4 * elem + 1], &mm->elems[4 * elem + 2], &mm->elems[4 * elem + 3]);
    if(cnt != 5) {
        printf("err: bad volumeelements\n");
        return 1;
    }
    return 0;
}
