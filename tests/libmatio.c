#include <stdio.h>
#include <unistd.h> /* for access() */
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

int main(int argc, char **argv)
{
    if( argc != 2 ) {
       fprintf(stderr, "missing .mat file\n");
       return 1;
    }
    mat_t * matfp = NULL;
    matvar_t * t = NULL;
    if( access( argv[1], F_OK ) == -1 ) {
       fprintf(stderr, "error: can't find %s\n", argv[1]);
       return 1;
    }
    printf("hello\n");
    matfp = Mat_Open(argv[1], MAT_ACC_RDONLY);
    if(matfp == NULL) {
       fprintf(stderr, "error: failed to open %s\n", argv[1]);
       return 2; /* bad fname */
    }
    int more_data = 1;
    int count = 0;
    while (more_data) {
        t = Mat_VarReadNextInfo(matfp);
        more_data = (t != NULL);
        if( more_data ) {
           // Mat_Rewind(matfp);
            matvar_t * tt = Mat_VarRead(matfp, t->name);
            Mat_VarPrint(tt, 1);
            count++;
        }
        Mat_VarFree(t);
    }
    if( count == 0 ) {
       fprintf(stderr, "no variables loaded from %s\n", argv[1]);
       return 3;
    }
    printf("goodbye\n");
    return 0;
}
