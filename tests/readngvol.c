#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf */
#include <libgen.h> /* basename */

#include "../src/readngvol.h"

int main(int argc, char ** argv)
{
    mesh m = {0};
    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol>\n", basename(argv[0]));
        return 1;
    }
    int ret = readngvol(argv[1], &m);
    if(ret != 0) {
        fprintf(stderr, "error: failed to load %s\n", argv[1]);
    }
    mesh_free(&m);
    return ret;
}
