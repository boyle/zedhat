#include <stdio.h>
#include <unistd.h> /* for access() */
#include <libgen.h> /* for basename() */
#include <string.h> /* for strcmp */
#include <matio.h>
#include "config.h"
/* testa.mat <- a
 * testb.mat <- b
 * testab.mat <- {a,b}
 * a =    1     2     3
 *        4     5     6
 * b = (1,1)        1
 *     (1,2)        2
 *     (2,2)        5 */
#define PROG basename(argv[0])
int main(int argc, char ** argv)
{
    if( argc != 3 ) {
        fprintf(stderr, "error: usage is '%s in.mat out.mat'\n", basename(argv[0]));
        return 1;
    }
    mat_t * in = NULL;
    mat_t * out = NULL;
    if( strcmp(argv[1], argv[2]) == 0 ) {
        fprintf(stderr, "%s error: in.mat and out.mat must differ: %s\n", PROG, argv[1]);
        return 1;
    }
    if( access( argv[1], F_OK ) == -1 ) {
        fprintf(stderr, "%s error: can't find %s\n", PROG, argv[1]);
        return 1;
    }
    if( access( argv[2], F_OK ) != -1 ) {
        if ( unlink(argv[2]) != 0 ) {
            fprintf(stderr, "%s error: %s already exists but refused to be deleted\n", PROG, argv[2]);
            return 1;
        }
    }
    printf("hello\n");
    in = Mat_Open(argv[1], MAT_ACC_RDONLY);
    if(in == NULL) {
        fprintf(stderr, "%s error: failed to open input %s\n", PROG, argv[1]);
        return 2; /* bad input filename */
    }
    out = Mat_Create(argv[2], "");
    if(out == NULL) {
        fprintf(stderr, "%s error: failed to open output %s\n", PROG, argv[2]);
        return 3; /* bad output filename */
    }
    int more_data = 1;
    int count = 0;
    while (more_data) {
        matvar_t * t = Mat_VarReadNextInfo(in);
        more_data = (t != NULL);
        if( more_data ) {
            t = Mat_VarRead(in, t->name);
            Mat_VarPrint(t, 1);
            Mat_VarWrite(out, t, t->compression);
            count++;
        }
        Mat_VarFree(t);
        t = NULL;
    }
    if( count == 0 ) {
        fprintf(stderr, "no variables loaded from %s\n", argv[1]);
        return 4; /* failed to load data */
    }
    Mat_Close(in);
    Mat_Close(out);
    printf("goodbye\n");
    return 0;
}
