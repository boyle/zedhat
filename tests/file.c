#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <string.h> /* strncmp, strlen */
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf */
#include <libgen.h> /* basename */
#include <zlib.h> /* cmocka: mocking zlib's gzopen, gzgets, gzclose, gzeof */

#include "cmocka.h"
#include "file.h"

int test_malloc_enabled;
void * __real__test_malloc(const size_t size, const char * file, const int line);
void * __wrap__test_malloc(size_t size)
{
    if(test_malloc_enabled && mock()) {
        return NULL;
    }
    else {
        return __real__test_malloc(size, __FILE__, __LINE__);
    }
}

#define MAGIC_FILE 0xdeadbeef
gzFile __real_gzopen(const char * path, const char * mode);
gzFile __wrap_gzopen(const char * path, const char * mode)
{
    if(strncmp(path, "test", 5) != 0) {
        return __real_gzopen(path, mode);
    }
    return (void *) mock();
}

char * __real_gzgets(gzFile file, char * buf, int len);
char * __wrap_gzgets(gzFile file, char * buf, int len)
{
    if((intptr_t) file != MAGIC_FILE) {
        return __real_gzgets(file, buf, len);
    }
    const char * str = mock_type(char *);
    strncpy(buf, str, len);
    return buf;
}

int __real_gzclose(gzFile file);
int __wrap_gzclose(gzFile file)
{
    if((intptr_t) file != MAGIC_FILE) {
        return __real_gzclose(file);
    }
    return mock();
}

int __real_gzeof(gzFile file);
int __wrap_gzeof(gzFile file)
{
    if((intptr_t) file != MAGIC_FILE) {
        return __real_gzeof(file);
    }
    return mock();
}

void mock_read(int n, int eof, int dim)
{
    if(n <= 0) {
        n = 16 + n;
    }
    will_return(__wrap_gzopen, MAGIC_FILE);
    will_return(__wrap_gzgets, "dimension\n");
    if(dim == 3) {
        will_return(__wrap_gzgets, "3\n");
    }
    else {
        will_return(__wrap_gzgets, "2\n");
    }
    will_return(__wrap_gzgets, "geomtype\n");
    will_return(__wrap_gzgets, "0\n");
    will_return(__wrap_gzgets, "volumeelements\n");
    will_return(__wrap_gzgets, "1\n");
    if(n >= 4) {
        if(dim == 3) {
            will_return(__wrap_gzgets, "1 4 4 2 6 8\n");
        }
        else {
            will_return(__wrap_gzgets, "1 3 4 2 6\n");
        }
        will_return(__wrap_gzgets, "points\n");
        will_return(__wrap_gzgets, "1\n");
    }
    if(n >= 5) {
        if(dim == 3) {
            will_return(__wrap_gzgets, "    0.0000000000000000      0.0000000000000000      0.0000000000000000\n");
        }
        else {
            will_return(__wrap_gzgets, "    0.0000000000000000      0.0000000000000000\n");
        }
        will_return(__wrap_gzgets, "surfaceelements\n");
        will_return(__wrap_gzgets, "1\n");
    }
    if(n >= 6) {
        if(dim == 3) {
            will_return(__wrap_gzgets, " 1 1 1 0 3 1 4 7\n");
        }
        else {
            will_return(__wrap_gzgets, " 1 1 1 0 2 1 4\n");
        }
    }
    will_return(__wrap_gzgets, "");
    will_return_count(__wrap_gzeof, 0, n);
    if(eof) {
        will_return(__wrap_gzeof, 1);
    }
    if(eof < 2) {
        will_return(__wrap_gzclose, Z_OK);
    }
    else {
        will_return(__wrap_gzclose, Z_STREAM_ERROR);
    }
}

void mock_ngvol_read(int n, int eof, int dim)
{
    will_return(__wrap_gzgets, "mesh3d\n");
    mock_read(n, eof, dim);
}

void mock_zh_read(int n, int eof, int dim)
{
    will_return(__wrap_gzgets, "zedhat\n");
    will_return(__wrap_gzgets, "format\n");
    will_return(__wrap_gzgets, "1\n");
    will_return(__wrap_gzeof, 0);
    /* ng.vol file, excepting optional first line 'mesh3d'*/
    mock_read(n, eof, dim);
}

static void test_happy (void ** state)
{
    int ret;
    (void) state; /* unused */
    model m = {{0}};
    int dim;
    for( dim = 2; dim <= 3; dim  += 1) {
        /* for a netgen file */
        printf("-- netgen %dd --\n", dim);
        mock_ngvol_read(6, 1, 3);
        will_return_count(__wrap__test_malloc, 0, 5);
        ret = readfile_ngvol("test", &m);
        assert_int_equal(ret, 0);
        model_free(&m);
        /* now repeat for a zedhat file */
        printf("-- zedhat %dd --\n", dim);
        mock_zh_read(6, 1, 3);
        will_return_count(__wrap__test_malloc, 0, 5);
        ret = readfile("test", &m);
        assert_int_equal(ret, 0);
        model_free(&m);
        /* now repeat for a zedhat file with optional fields: stimmeas */
        printf("-- zedhat %dd (optional fields) --\n", dim);
        will_return_count(__wrap_gzeof, 0, 4);
        mock_zh_read(6, 1, 3);
        will_return(__wrap_gzgets, "stimmeas\n");
        will_return(__wrap_gzgets, "1\n");
        will_return(__wrap_gzgets, "1 2 3 4\n");
        will_return(__wrap_gzgets, "data\n");
        will_return(__wrap_gzgets, "1 2\n");
        will_return(__wrap_gzgets, " 1.1 -2.2\n");
        will_return(__wrap_gzgets, "parameters\n");
        will_return(__wrap_gzgets, "1 3\n");
        will_return(__wrap_gzgets, "-1.1 2.2 +3.3 4.4e-2\n");
        will_return(__wrap_gzgets, "hyperparameter\n");
        will_return(__wrap_gzgets, " 1.1\n");
        will_return_count(__wrap__test_malloc, 0, 5 + 3);
        ret = readfile("test", &m);
        assert_int_equal(ret, 0);
        model_free(&m);
    }
}

static void test_mangle_fail1 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return_count(__wrap_gzeof, 0, 1);
    mock_zh_read(6, 0, 3);
    will_return(__wrap_gzgets, "hyperparameter\n");
    will_return(__wrap_gzgets, " asdf\n");
    will_return_count(__wrap__test_malloc, 0, 5);
    ret = readfile("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail1 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap__test_malloc, 1);
    will_return(__wrap__test_malloc, 0);
    mock_ngvol_read(3, 0, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail2 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 1);
    mock_ngvol_read(3, 0, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail3 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 1);
    mock_ngvol_read(4, 0, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail4 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 1);
    will_return(__wrap__test_malloc, 0);
    mock_ngvol_read(5, 0, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail5 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 0);
    will_return(__wrap__test_malloc, 1);
    mock_ngvol_read(5, 0, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail6 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap_gzeof, 0);
    mock_zh_read(6, 0, 3);
    will_return(__wrap_gzgets, "stimmeas\n");
    will_return(__wrap_gzgets, "1\n");
    will_return(__wrap_gzgets, "\n");
    will_return_count(__wrap__test_malloc, 0, 5);
    will_return(__wrap__test_malloc, 1);
    ret = readfile("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail7 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap_gzeof, 0);
    mock_zh_read(6, 0, 3);
    will_return(__wrap_gzgets, "data\n");
    will_return(__wrap_gzgets, "2 6\n");
    will_return(__wrap_gzgets, "\n");
    will_return_count(__wrap__test_malloc, 0, 5);
    will_return(__wrap__test_malloc, 1);
    ret = readfile("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_malloc_fail8 (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return(__wrap_gzeof, 0);
    mock_zh_read(6, 0, 3);
    will_return(__wrap_gzgets, "parameters\n");
    will_return(__wrap_gzgets, "2 6\n");
    will_return(__wrap_gzgets, "\n");
    will_return_count(__wrap__test_malloc, 0, 5);
    will_return(__wrap__test_malloc, 1);
    ret = readfile("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_gzclose_fail (void ** state)
{
    (void) state; /* unused */
    int ret;
    model m = {{0}};
    will_return_count(__wrap__test_malloc, 0, 5);
    mock_ngvol_read(6, 2, 3);
    ret = readfile_ngvol("test", &m);
    assert_int_equal(ret, 1);
    model_free(&m);
}

static void test_null_fail (void ** state)
{
    (void) state; /* unused */
    assert_int_equal(readfile_ngvol("test", NULL), 1);
    assert_int_equal(readfile("test", NULL), 1);
}

int main(int argc, char ** argv)
{
    test_malloc_enabled = 0;
    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol|test>\n", basename(argv[0]));
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
            cmocka_unit_test(test_malloc_fail6),
            cmocka_unit_test(test_malloc_fail7),
            cmocka_unit_test(test_malloc_fail8),
            cmocka_unit_test(test_gzclose_fail),
            cmocka_unit_test(test_null_fail),
        };
        return cmocka_run_group_tests(tests, NULL, NULL);
    }
    model m = {{0}};
    int len = strlen(filename);
    int ret;
    if(((len > 4) && (strcmp(&(filename[len - 4]), ".vol") == 0)) ||
       ((len > 7) && (strcmp(&(filename[len - 7]), ".vol.gz") == 0))) {
        ret = readfile_ngvol(filename, &m);
    }
    else {
        ret = readfile(filename, &m);
    }
    if(ret != 0) {
        fprintf(stderr, "error: failed to load %s\n", argv[1]);
    }
    model_free(&m);
    return ret;
}
