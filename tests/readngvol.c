#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf */
#include <libgen.h> /* basename */

#include "../src/readngvol.h"

int main(int argc, char ** argv)
{
    struct mesh m = {0};

    if(argc != 2) {
        fprintf(stderr, "usage: %s <netgen.vol>\n", basename(argv[0]));
        return 1;
    }

    int ret = readngvol(argv[1], &m);
    if(ret != 0)
       fprintf(stderr,"error: failed to load %s\n",argv[1]);

    free(m.nodes); m.nodes = NULL;
    free(m.elems); m.elems = NULL;
    free(m.matidx); m.matidx = NULL;
    free(m.surfaceelems); m.surfaceelems = NULL;
    free(m.bc); m.bc = NULL;

    return ret;
}
