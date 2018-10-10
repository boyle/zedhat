#include <string.h> /* memset */
#include <stdlib.h> /* free */

#include "model.h"

void mesh_init(struct mesh * m)
{
    memset(m, 0, sizeof(struct mesh));
}

void mesh_free(struct mesh * m)
{
    free(m->nodes);
    free(m->elems);
    free(m->matidx);
    free(m->surfaceelems);
    free(m->bc);
    mesh_init(m);
}
