/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h> /* printf */
#include <stddef.h> /* offsetof */
#include <stdint.h> /* SIZE_MAX */
#include <string.h> /* strlen, strncat */
#include <ctype.h> /* isspace */
#include <zlib.h> /* gzopen, gzgets, gzeof, gzclose */
#include <assert.h> /* assert */
#include <limits.h> /* MAX_INT */

#include "model.h"
#include "file.h"

#ifdef UNIT_TESTING
extern gzFile _mock_gzopen(const char * path, const char * mode);
extern char * _mock_gzgets(gzFile file, char * buf, int len);
extern int _mock_gzclose(gzFile file);
extern int _mock_gzeof(gzFile file);
#define gzopen(path, mode) _mock_gzopen(path, mode)
#define gzgets(file, buf, len) _mock_gzgets(file, buf, len)
#define gzclose(file) _mock_gzclose(file)
#define gzeof(file) _mock_gzeof(file)
#endif

#define SUCCESS 1
#define FAILURE 0

#define MAXCHAR 1024
static void gzreadnext(gzFile F, char data[], int n)
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
enum fileformat_field {EOL, TEST, FMT, MESH, HP, STIM, DATA, PARAM, DIM, TYPE, SURF, NODE, ELEM, MAP, ZC};
enum fileformat_req {OPTIONAL, REQUIRED, REQUIRED_FIRST};
typedef struct {
    enum fileformat_field field;
    char section [LEN];
    int (*sscanf_funcptr)(const char *, model *, mesh *, int i);
    enum fileformat_req req;
    int cnt;
} fileformat;

static int check_req(const fileformat * f, enum fileformat_req needs)
{
    for(int i = 0; f[i].section[0]; i++) {
        if((f[i].req >= needs) && (!f[i].cnt)) {
            return i + 1;
        }
    }
    return 0;
}

static int ooo_section(const char * section)
{
    printf("err: expect dimension before %s\n", section);
    return -1;
}

static int bad_section(const char * section)
{
    printf("err: bad %s\n", section);
    return -1;
}

static int long_section(const char * section)
{
    printf("err: %s too long\n", section);
    return -1;
}

static int bad_malloc(const char * section)
{
    printf("err: out of memory %s\n", section);
    return -1;
}

static int sscanf_size(const char * data, const char * section, enum fileformat_field field, model * m, mesh ** mm)
{
    assert(m != NULL);
    assert(mm != NULL);
    assert(*mm != NULL);
    static const char fwd[] = "forward\n";
    static const char rec[] = "reconstruction\n";
    int cnt = 0;
    int expect = 1;
    int n[2] = {1, 1};
    double d = 0;
    switch(field) {
    case MESH:
        expect = 0;
        break;
    case HP:
        expect = -1; /* double */
        break;
    case PARAM:
    case DATA:
        expect = 2; /* 2 ints */
        break;
    default:
        expect = 1; /* 1 int */
        break;
    }
    if(expect < 0) {
        expect = -expect;
        cnt = sscanf(data, " %lf\n", &d);
        if(cnt == expect) {
            printf("%s %lg\n", section, d);
        }
    }
    else if(expect > 0) {
        for(int j = 0; j < expect; j++) {
            int inc = 0;
            cnt += sscanf(data, " %d%n", &(n[j]), &inc);
            data += inc;
        }
        if(cnt == expect) {
            if(n[1] > 1) {
                printf("%s %dx%d\n", section, n[0], n[1]);
            }
            else {
                printf("%s %d\n", section, n[0]);
            }
        }
    }
    const int min_rows = (field == TYPE) ? 0 : 1;
    if((cnt != expect) || (n[0] < min_rows) || (n[1] < 1) || (d < 0.0)) {
        return bad_section(section);
    }
    int rows = n[0];
    if(rows > INT_MAX / (sizeof(double) > sizeof(int) ? sizeof(double) : sizeof(int)) / 10) {
        return long_section(section);
    }
    switch(field) {
    case FMT:
        if(n[0] != 1) {
            return bad_section(section);
        }
        rows = 0;
        break;
    case MESH:
        if(strncmp(data, fwd, strlen(fwd)) == 0) {
            *mm = &(m->fwd);
            printf("%s %s", section, fwd);
        }
        else if(strncmp(data, rec, strlen(rec)) == 0) {
            *mm = &(m->rec);
            printf("%s %s", section, rec);
        }
        else {
            return bad_section(section);
        }
        rows = 0;
        break;
    case STIM:
        if(!set_model_stimmeas(m, n[0])) {
            return bad_malloc(section);
        }
        break;
    case DATA:
        if(!set_model_data(m, n[0], n[1])) {
            return bad_malloc(section);
        }
        break;
    case PARAM:
        if(!set_model_params(m, n[0], n[1])) {
            return bad_malloc(section);
        }
        break;
    case HP:
        rows = 1;
        break;
    case DIM:
        if(n[0] != 2 && n[0] != 3) {
            return bad_section(section);
        }
        set_mesh_dim(*mm, n[0]);
        rows = 0;
        break;
    case TYPE:
        if(n[0] != 0) {
            return bad_section(section);
        }
        rows = 0;
        break;
    case SURF:
        if((*mm)->dim == 0) {
            return ooo_section(section);
        }
        if(!set_mesh_surfaceelems(*mm, n[0])) {
            return bad_malloc(section);
        }
        break;
    case NODE:
        if((*mm)->dim == 0) {
            return ooo_section(section);
        }
        if(!set_mesh_nodes(*mm, n[0])) {
            return bad_malloc(section);
        }
        break;
    case ELEM:
        if((*mm)->dim == 0) {
            return ooo_section(section);
        }
        if(!set_mesh_elems(*mm, n[0])) {
            return bad_malloc(section);
        }
        break;
    case MAP:
        if(!set_mesh_pmap(*mm, n[0])) {
            return bad_malloc(section);
        }
        break;
    case ZC:
        if(!set_model_zc(m, n[0])) {
            return bad_malloc(section);
        }
        break;
    default:  /* LCOV_EXCL_LINE */
        assert(0);  /* LCOV_EXCL_LINE */
    }
    return rows;
}


int readfile_loop(const char filename[], model * m, fileformat * f_list)
{
    if(m == NULL) {
        return FAILURE;
    }
    reset_model(m);
    int ret = FAILURE;
    char data[MAXCHAR];
    gzFile F = gzopen(filename, "r");
    if(!F) {
        printf("err: failed to open %s\n", filename);
        goto __quit;
    }
    printf("reading %s\n", filename);
    mesh * mm = &(m->fwd); /* default: forward mesh */
    int first = 1;
    do {
        gzreadnext(F, data, MAXCHAR);
        fileformat * f = &(f_list[0]);
        while(f->field != EOL) {
            char section_n [LEN + 1];
            strncpy(section_n, f->section, LEN + 1);
            section_n[LEN - 1] = '\0';
            const int eol = strlen(section_n);
            strcpy(&section_n[eol], "\n\0");
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
            if(f->field == TEST) {
                printf("%s\n", f->section);
                break;
            }
            gzreadnext(F, data, MAXCHAR);
            int n = sscanf_size(data, f->section, f->field, m, &mm);
            if(n < 0) {
                goto __quit;
            }
            for(int i = 0; i < n; i += 1) {
                gzreadnext(F, data, MAXCHAR);
                int rets = (*f->sscanf_funcptr)(data, m, mm, i);
                if( rets == FAILURE ) {
                    bad_section(f->section);
                    goto __quit;
                }
            }
        }
    } while(!gzeof(F));
    int chk = check_req(f_list, REQUIRED);
    if(chk) {
        printf("err: %s: missing %s\n", filename, (f_list[chk - 1]).section);
        goto __quit;
    }
    ret = SUCCESS;
__quit:
    if(F) {
        if(gzclose(F) != Z_OK) {
            printf("err: failed to close %s\n", filename);
            ret = FAILURE;
        }
    }
    return ret;
}

static int readngvol_surfaceelements(const char * data, model * m, mesh * mm, const int se)
{
    assert(mm != NULL);
    assert(mm->surfaceelems != NULL);
    assert(mm->bc != NULL);
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, "%*d %d %*d %*d 3 %d %d %d\n", &mm->bc[se], &mm->surfaceelems[3 * se + 0], &mm->surfaceelems[3 * se + 1], &mm->surfaceelems[3 * se + 2]);
    }
    else {
        cnt = sscanf(data, "%*d %d %*d %*d 2 %d %d\n", &mm->bc[se], &mm->surfaceelems[2 * se + 0], &mm->surfaceelems[2 * se + 1]);
    }
    return (cnt == mm->dim + 1) ? SUCCESS : FAILURE;
}

static int readngvol_points(const char * data, model * m, mesh * mm, const int node)
{
    assert(mm != NULL);
    assert(mm->nodes != NULL);
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, " %lf %lf %lf\n", &mm->nodes[3 * node + 0], &mm->nodes[3 * node + 1], &mm->nodes[3 * node + 2]);
    }
    else {
        cnt = sscanf(data, " %lf %lf\n", &mm->nodes[2 * node + 0], &mm->nodes[2 * node + 1]);
    }
    return (cnt == mm->dim) ? SUCCESS : FAILURE;
}

static int readngvol_volumeelements(const char * data, model * m, mesh * mm, const int elem)
{
    assert(mm != NULL);
    assert(mm->elems != NULL);
    assert(mm->matidx != NULL);
    int cnt;
    if(mm->dim == 3) {
        cnt = sscanf(data, "%d 4 %d %d %d %d\n", &mm->matidx[elem], &mm->elems[4 * elem + 0], &mm->elems[4 * elem + 1], &mm->elems[4 * elem + 2], &mm->elems[4 * elem + 3]);
    }
    else {
        cnt = sscanf(data, "%d 3 %d %d %d\n", &mm->matidx[elem], &mm->elems[3 * elem + 0], &mm->elems[3 * elem + 1], &mm->elems[3 * elem + 2]);
    }
    return (cnt == mm->dim + 2) ? SUCCESS : FAILURE;
}

static int readzh_stimmeas(const char * data, model * m, mesh * mm, const int i)
{
    assert(m != NULL);
    assert(m->stimmeas != NULL);
    int cnt = sscanf(data, "%d %d %d %d %lf\n",
                     &m->stimmeas[4 * i + 0], &m->stimmeas[4 * i + 1],
                     &m->stimmeas[4 * i + 2], &m->stimmeas[4 * i + 3],
                     &m->measgain[i]);
    return (cnt == 5) ? SUCCESS : FAILURE;
}

static int readzh_data(const char * data, model * m, mesh * mm, const int i)
{
    assert(m != NULL);
    assert(m->data != NULL);
    const int rows = m->n_data[0];
    const int cols = m->n_data[1];
    int cnt = 0;
    for(int j = 0; j < cols; j++) {
        int inc = 0;
        cnt += sscanf(data, " %lf %n", &(m->data[i + rows * j]), &inc);
        data += inc;
    }
    return (cnt == cols) ? SUCCESS : FAILURE;
}

static int readzh_parametermap(const char * data, model * m, mesh * mm, const int i)
{
    assert(mm != NULL);
    int cnt = sscanf(data, " %d %d %lf\n", &(mm->pmap_elem[i]), &(mm->pmap_param[i]), &(mm->pmap_frac[i]));
    return (cnt == 3) ? SUCCESS : FAILURE;
}

static int readzh_params(const char * data, model * m, mesh * mm, const int i)
{
    assert(m != NULL);
    assert(m->params != NULL);
    const int rows = m->n_params[0];
    const int cols = m->n_params[1];
    int cnt = 0;
    for(int j = 0; j < cols; j++) {
        int inc = 0;
        cnt += sscanf(data, " %lf%n", &(m->params[i + rows * j]), &inc);
        data += inc;
    }
    return (cnt == cols) ? SUCCESS : FAILURE;
}


static int readzh_hp(const char * data, model * m, mesh * mm, const int i)
{
    assert(m != NULL);
    int cnt = sscanf(data, " %lf\n", &(m->hp));
    return (cnt == 1) ? SUCCESS : FAILURE;
}

static int readzh_zc(const char * data, model * m, mesh * mm, const int i)
{
    assert(mm != NULL);
    int cnt = sscanf(data, " %d %lf\n", &(m->zcbc[i]), &(m->zc[i]));
    return (cnt == 2) ? SUCCESS : FAILURE;
}


int readfile_ngvol(const char filename[], model * m)
{
    fileformat ngvol [] = {
        {TEST, "mesh3d", NULL, REQUIRED_FIRST},
        {DIM,  "dimension", NULL, REQUIRED},
        {TYPE, "geomtype", NULL, REQUIRED},
        /* GEOM_STL=11: surfaceelementsgi
        * GEOM_OCC=12, GEOM_ACIS=13: surfaceelementsuv
        * default: surfaceelements
        */
        {SURF, "surfaceelements", &readngvol_surfaceelements, REQUIRED},
        {NODE, "points", &readngvol_points, REQUIRED},
        {ELEM, "volumeelements", &readngvol_volumeelements, REQUIRED},
        { 0 }
    };
    return readfile_loop(filename, m, ngvol);
}

int readfile(const char filename[], model * m)
{
    fileformat zh_format1 [] = {
        {TEST, "zedhat", NULL, REQUIRED_FIRST},
        {FMT, "format", NULL, REQUIRED},
        {MESH, "modeltype", NULL, OPTIONAL},
        {STIM, "stimmeas", &readzh_stimmeas, OPTIONAL},
        {DATA, "data", &readzh_data, OPTIONAL},
        {PARAM, "parameters", &readzh_params, OPTIONAL},
        {HP,   "hyperparameter", &readzh_hp, OPTIONAL},
        /* netgen .vol format */
        {DIM,  "dimension", NULL, REQUIRED},
        {TYPE, "geomtype", NULL, REQUIRED},
        {SURF, "surfaceelements", &readngvol_surfaceelements, REQUIRED},
        {NODE, "points", &readngvol_points, REQUIRED},
        {ELEM, "volumeelements", &readngvol_volumeelements, REQUIRED},
        {MAP, "parametermap", &readzh_parametermap, OPTIONAL},
        {ZC, "contactimpedances", &readzh_zc, OPTIONAL},
        { 0 }
    };
    if(!readfile_loop(filename, m, zh_format1)) {
        return FAILURE;
    }
    return SUCCESS;
}
