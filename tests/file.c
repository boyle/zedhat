/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <string.h> /* strncmp, strlen */
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf, sprintf */
#include <limits.h> /* INT_MAX */
#include <libgen.h> /* basename */
#include <zlib.h> /* cmocka: mocking zlib's gzopen, gzgets, gzclose, gzeof */

#include "file.h"
#include "cmocka.h"

int test_malloc_enabled;
void * _mock_test_malloc(const size_t size, const char * file, const int line)
{
    if(test_malloc_enabled && mock()) {
        printf("  malloc x\n");
        return NULL;
    }
    else {
        if(test_malloc_enabled) {
            printf("  malloc ✓\n");
        }
        return _test_malloc(size, file, line);
    }
}

#define MAGIC_FILE 0xdeadbeef
gzFile _mock_gzopen(const char * path, const char * mode)
{
    if(strncmp(path, "unittest", 5) != 0) {
        return gzopen(path, mode);
    }
    return (void *) mock();
}

char * _mock_gzgets(gzFile file, char * buf, int len)
{
    if((intptr_t) file != MAGIC_FILE) {
        return gzgets(file, buf, len);
    }
    const char * str = mock_type(char *);
    strncpy(buf, str, len);
    return buf;
}

int _mock_gzclose(gzFile file)
{
    if((intptr_t) file != MAGIC_FILE) {
        return gzclose(file);
    }
    return mock();
}

int _mock_gzeof(gzFile file)
{
    if((intptr_t) file != MAGIC_FILE) {
        return gzeof(file);
    }
    return mock();
}

void mock_read(int n, int eof, int dim, int pmap)
{
    if(n <= 0) {
        n = 16 + n;
    }
    will_return(_mock_gzopen, MAGIC_FILE);
    will_return(_mock_gzgets, "dimension\n");
    if(dim == 3) {
        will_return(_mock_gzgets, "3\n");
    }
    else {
        will_return(_mock_gzgets, "2\n");
    }
    will_return(_mock_gzgets, "geomtype\n");
    will_return(_mock_gzgets, "0\n");
    will_return(_mock_gzgets, "volumeelements\n");
    will_return(_mock_gzgets, "1\n");
    if(n >= 4) {
        if(dim == 3) {
            will_return(_mock_gzgets, "1 4 4 2 6 8\n");
        }
        else {
            will_return(_mock_gzgets, "1 3 4 2 6\n");
        }
        will_return(_mock_gzgets, "points\n");
        will_return(_mock_gzgets, "1\n");
    }
    if(n >= 5) {
        if(dim == 3) {
            will_return(_mock_gzgets, "    0.0000000000000000      0.0000000000000000      0.0000000000000000\n");
        }
        else {
            will_return(_mock_gzgets, "    0.0000000000000000      0.0000000000000000\n");
        }
        will_return(_mock_gzgets, "surfaceelements\n");
        will_return(_mock_gzgets, "1\n");
    }
    if(n >= 6) {
        if(dim == 3) {
            will_return(_mock_gzgets, " 1 1 1 0 3 1 4 7\n");
        }
        else {
            will_return(_mock_gzgets, " 1 1 1 0 2 1 4\n");
        }
        if(pmap) {
            will_return(_mock_gzgets, "parametermap\n");
            will_return(_mock_gzgets, "1\n");
        }
    }
    if(n >= 7) {
        will_return(_mock_gzgets, " 1 1 1.01\n");
    }
    if(eof) {
        will_return(_mock_gzgets, "");
    }
    will_return_count(_mock_gzeof, 0, n);
    if(eof) {
        will_return(_mock_gzeof, 1);
    }
    if(eof < 2) {
        will_return(_mock_gzclose, Z_OK);
    }
    else {
        will_return(_mock_gzclose, Z_STREAM_ERROR);
    }
}

void mock_ngvol_read(int n, int eof, int dim)
{
    assert_int_equal(n <= 6, 1);
    will_return(_mock_gzgets, "mesh3d\n");
    mock_read(n, eof, dim, 0);
}

void mock_zh_read(int n, int eof, int dim)
{
    assert_int_equal(n <= 7, 1); /* TODO only loading fwd models so far, no rec models */
    will_return(_mock_gzgets, "zedhat\n");
    will_return(_mock_gzgets, "format\n");
    will_return(_mock_gzgets, "1\n");
    will_return(_mock_gzeof, 0);
    will_return(_mock_gzgets, "modeltype\n");
    will_return(_mock_gzgets, "forward\n");
    will_return(_mock_gzeof, 0);
    /* ng.vol file, excepting optional first line 'mesh3d'*/
    mock_read(n, (n <= 7) ? eof : 0, dim, 1);
    if(n > 9) {
        will_return(_mock_gzgets, "modeltype\n");
        will_return(_mock_gzgets, "reconstruction\n");
        will_return(_mock_gzeof, 0);
        /* ng.vol file, excepting optional first line 'mesh3d'*/
        mock_read(n - 9, eof, dim, 1);
    }
}

static void test_happy (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    for(int dim = 2; dim <= 3; dim  += 1) {
        /* for a netgen file */
        printf("-- netgen %dd --\n", dim);
        will_return_count(_mock_test_malloc, 0, 5);
        mock_ngvol_read(6, 1, 3);
        ret = readfile_ngvol("unittest", m);
        assert_int_equal(ret, 1);
        /* now repeat for a zedhat file */
        printf("-- zedhat %dd --\n", dim);
        will_return_count(_mock_test_malloc, 0, 5 + 3);
        mock_zh_read(7, 1, 3);
        ret = readfile("unittest", m);
        assert_int_equal(ret, 1);
        /* now repeat for a zedhat file with optional fields: stimmeas */
        printf("-- zedhat %dd (optional fields) --\n", dim);
        will_return_count(_mock_test_malloc, 0, 5 + 3 + 3 + 3);
        will_return_count(_mock_gzeof, 0, 5);
        mock_zh_read(7, 1, 3);
        will_return(_mock_gzgets, "stimmeas\n");
        will_return(_mock_gzgets, "1\n");
        will_return(_mock_gzgets, "1 2 3 4 2.0\n");
        will_return(_mock_gzgets, "data\n");
        will_return(_mock_gzgets, "1 2\n");
        will_return(_mock_gzgets, " 1.1 -2.2\n");
        will_return(_mock_gzgets, "parameters\n");
        will_return(_mock_gzgets, "1 3\n");
        will_return(_mock_gzgets, "-1.1 2.2 +3.3 4.4e-2\n");
        will_return(_mock_gzgets, "hyperparameter\n");
        will_return(_mock_gzgets, " 1\n");
        will_return(_mock_gzgets, " 1.1\n");
        will_return(_mock_gzgets, "contactimpedances\n");
        will_return(_mock_gzgets, " 1\n");
        will_return(_mock_gzgets, "1 10.0\n");
        ret = readfile("unittest", m);
        assert_int_equal(ret, 1);
    }
    m = free_model(m);
}

static void test_mangle_fail_hp_mangled (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(7, 0, 3);
    will_return(_mock_gzgets, "hyperparameter\n");
    will_return(_mock_gzgets, " asdf\n");
    will_return_count(_mock_test_malloc, 0, 8);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_mangle_modeltype (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(7, 0, 3);
    will_return(_mock_gzgets, "modeltype\n");
    will_return(_mock_gzgets, "asdf\n");
    will_return_count(_mock_test_malloc, 0, 8);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_ngvol_malloc_fails (void ** state)
{
    for(int i = 0; i < 5; i++) {
        int n = 0;
        switch(i) {
        case 0:
        case 1:
            n = 3;
            break;
        case 2:
            n = 4;
            break;
        case 3:
        case 4:
            n = 5;
            break;
        default: assert_int_equal(1, 0);
        }
        printf("i=%d, n=%d\n", i, n);
        will_return(_mock_test_malloc, 0);
        model * m = malloc_model();
        if(i > 0) {
            will_return_count(_mock_test_malloc, 0, i);
        }
        will_return(_mock_test_malloc, 1);
        mock_ngvol_read(n, 0, 3);
        int ret = readfile_ngvol("unittest", m);
        assert_int_equal(ret, 0);
        free_model(m);
    }
}

static void test_zh_malloc_fails (void ** state)
{
    for(int i = 0; i < 8; i++) {
        int n = 0;
        switch(i) {
        case 0:
        case 1:
            n = 3;
            break;
        case 2:
            n = 4;
            break;
        case 3:
        case 4:
            n = 5;
            break;
        case 5:
        case 6:
        case 7:
            n = 6;
            break;
        default: assert_int_equal(1, 0);
        }
        printf("i=%d, n=%d\n", i, n);
        will_return(_mock_test_malloc, 0);
        model * m = malloc_model();
        if(i > 0) {
            will_return_count(_mock_test_malloc, 0, i);
        }
        will_return(_mock_test_malloc, 1);
        mock_zh_read(n, 0, 3);
        int ret = readfile("unittest", m);
        assert_int_equal(ret, 0);
        free_model(m);
    }
}

static void test_zh_malloc_fail_optional (void ** state)
{
    for(int i = 0; i < 6; i++) {
        printf("i=%d\n", i);
        int ret;
        will_return(_mock_test_malloc, 0);
        model * m = malloc_model();
        mock_zh_read(7, 0, 3);
        int extra_malloc = 0;
        switch(i) {
        case 1:
            extra_malloc = 1;
        case 0:
            will_return(_mock_gzgets, "stimmeas\n");
            will_return(_mock_gzgets, "1\n");
            break;
        case 2:
            will_return(_mock_gzgets, "data\n");
            will_return(_mock_gzgets, "2 6\n");
            break;
        case 3:
            will_return(_mock_gzgets, "parameters\n");
            will_return(_mock_gzgets, "2 6\n");
            break;
        case 5:
            extra_malloc = 1;
        case 4:
            will_return(_mock_gzgets, "contactimpedances\n");
            will_return(_mock_gzgets, "1\n");
            break;
        default: fail();
        }
        will_return_count(_mock_test_malloc, 0, 8 + extra_malloc);
        will_return(_mock_test_malloc, 1);
        ret = readfile("unittest", m);
        assert_int_equal(ret, 0);
        free_model(m);
    }
}


static void test_long_fail (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(7, 0, 3);
    will_return(_mock_gzgets, "parameters\n");
    char too_many[1024];
    sprintf(too_many, "%d 6\n", INT_MAX - 1);
    will_return(_mock_gzgets, too_many);
    will_return_count(_mock_test_malloc, 0, 8);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_gzclose_fail (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return_count(_mock_test_malloc, 0, 5);
    mock_ngvol_read(6, 2, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_null_fail (void ** state)
{
    (void) state; /* unused */
    assert_int_equal(readfile_ngvol("unittest", NULL), 0);
    assert_int_equal(readfile("unittest", NULL), 0);
}

int main(int argc, char ** argv)
{
    test_malloc_enabled = 0;
    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol|unittest>\n", basename(argv[0]));
        return 1;
    }
    char * filename = argv[1];
    if(strncmp(filename, "unittest", 9) == 0) {
        test_malloc_enabled = 1;
        const struct CMUnitTest tests[] = {
            cmocka_unit_test(test_happy),
            cmocka_unit_test(test_mangle_fail_hp_mangled),
            cmocka_unit_test(test_mangle_modeltype),
            cmocka_unit_test(test_ngvol_malloc_fails),
            cmocka_unit_test(test_zh_malloc_fails),
            cmocka_unit_test(test_zh_malloc_fail_optional),
            cmocka_unit_test(test_long_fail),
            cmocka_unit_test(test_gzclose_fail),
            cmocka_unit_test(test_null_fail),
        };
        return cmocka_run_group_tests(tests, NULL, NULL);
    }
    model * m = malloc_model();
    int len = strlen(filename);
    int ret;
    if(((len > 4) && (strcmp(&(filename[len - 4]), ".vol") == 0)) ||
       ((len > 7) && (strcmp(&(filename[len - 7]), ".vol.gz") == 0))) {
        ret = !readfile_ngvol(filename, m);
    }
    else {
        ret = !readfile(filename, m);
    }
    if(ret) {
        fprintf(stderr, "error: failed to load %s\n", argv[1]);
    }
    free_model(m);
    return ret;
}
