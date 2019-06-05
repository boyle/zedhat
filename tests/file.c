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

void mock_read(int n, int eof, int dim)
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
    will_return(_mock_gzgets, "mesh3d\n");
    mock_read(n, eof, dim);
}

void mock_zh_read(int n, int eof, int dim)
{
    will_return(_mock_gzgets, "zedhat\n");
    will_return(_mock_gzgets, "format\n");
    will_return(_mock_gzgets, "1\n");
    will_return(_mock_gzeof, 0);
    /* ng.vol file, excepting optional first line 'mesh3d'*/
    mock_read(n, eof, dim);
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
        will_return_count(_mock_test_malloc, 0, 5);
        mock_zh_read(6, 1, 3);
        ret = readfile("unittest", m);
        assert_int_equal(ret, 1);
        /* now repeat for a zedhat file with optional fields: stimmeas */
        printf("-- zedhat %dd (optional fields) --\n", dim);
        will_return_count(_mock_test_malloc, 0, 5 + 4);
        will_return_count(_mock_gzeof, 0, 4);
        mock_zh_read(6, 1, 3);
        will_return(_mock_gzgets, "stimmeas\n");
        will_return(_mock_gzgets, "1\n");
        will_return(_mock_gzgets, "1 2 3 4\n");
        will_return(_mock_gzgets, "data\n");
        will_return(_mock_gzgets, "1 2\n");
        will_return(_mock_gzgets, " 1.1 -2.2\n");
        will_return(_mock_gzgets, "parameters\n");
        will_return(_mock_gzgets, "1 3\n");
        will_return(_mock_gzgets, "-1.1 2.2 +3.3 4.4e-2\n");
        will_return(_mock_gzgets, "hyperparameter\n");
        will_return(_mock_gzgets, " 1.1\n");
        ret = readfile("unittest", m);
        assert_int_equal(ret, 1);
    }
    m = free_model(m);
}

static void test_mangle_fail1 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(6, 0, 3);
    will_return(_mock_gzgets, "hyperparameter\n");
    will_return(_mock_gzgets, " asdf\n");
    will_return_count(_mock_test_malloc, 0, 5);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail1 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_test_malloc, 1);
    will_return(_mock_test_malloc, 0);
    mock_ngvol_read(3, 0, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail2 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 1);
    mock_ngvol_read(3, 0, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail3 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 1);
    mock_ngvol_read(4, 0, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail4 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 1);
    will_return(_mock_test_malloc, 0);
    mock_ngvol_read(5, 0, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail5 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 1);
    mock_ngvol_read(5, 0, 3);
    ret = readfile_ngvol("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail6a (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(6, 0, 3);
    will_return(_mock_gzgets, "stimmeas\n");
    will_return(_mock_gzgets, "1\n");
    will_return_count(_mock_test_malloc, 0, 5);
    will_return(_mock_test_malloc, 1);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail6b (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    will_return(_mock_gzeof, 0);
    mock_zh_read(6, 1, 3);
    will_return(_mock_gzgets, "stimmeas\n");
    will_return(_mock_gzgets, "1\n");
    will_return(_mock_gzgets, "1 2 3 4\n");
    will_return_count(_mock_test_malloc, 0, 6);
    will_return(_mock_test_malloc, 1);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail7 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(6, 0, 3);
    will_return(_mock_gzgets, "data\n");
    will_return(_mock_gzgets, "2 6\n");
    will_return_count(_mock_test_malloc, 0, 5);
    will_return(_mock_test_malloc, 1);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_malloc_fail8 (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(6, 0, 3);
    will_return(_mock_gzgets, "parameters\n");
    will_return(_mock_gzgets, "2 6\n");
    will_return_count(_mock_test_malloc, 0, 5);
    will_return(_mock_test_malloc, 1);
    ret = readfile("unittest", m);
    assert_int_equal(ret, 0);
    free_model(m);
}

static void test_long_fail (void ** state)
{
    int ret;
    will_return(_mock_test_malloc, 0);
    model * m = malloc_model();
    mock_zh_read(6, 0, 3);
    will_return(_mock_gzgets, "parameters\n");
    char too_many[1024];
    sprintf(too_many, "%d 6\n", INT_MAX - 1);
    will_return(_mock_gzgets, too_many);
    will_return_count(_mock_test_malloc, 0, 5);
    will_return(_mock_test_malloc, 1);
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
            cmocka_unit_test(test_mangle_fail1),
            cmocka_unit_test(test_malloc_fail1),
            cmocka_unit_test(test_malloc_fail2),
            cmocka_unit_test(test_malloc_fail3),
            cmocka_unit_test(test_malloc_fail4),
            cmocka_unit_test(test_malloc_fail5),
            cmocka_unit_test(test_malloc_fail6a),
            cmocka_unit_test(test_malloc_fail6b),
            cmocka_unit_test(test_malloc_fail7),
            cmocka_unit_test(test_malloc_fail8),
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
