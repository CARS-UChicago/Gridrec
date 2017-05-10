#include <stdio.h>
#include <stdlib.h>

/* Memory allocation routines */
/*
 * The calloc routines assume all pointers are of the same size as
 * a char pointer. In order to guarantee that the matrices occupy
 * contiguous storage,  the higher dimension allocators do not call
 * lower dimension ones repeatedly. 
 */

void alloc_error_print(int n1, int n2)
{
    fprintf(stderr, "Allocation error in callocD routines\n");
    fprintf(stderr, "trying to allocate %d elements of size %d\n", n1, n2);
}

void **
calloc2D(int nr, int nc, int size)
{
    char **p;
    int i;
    
    p = (char **)calloc(nr, sizeof(char *));
    if ( !p ) {
	alloc_error_print(nr, sizeof(char *));
	return NULL;
    }
    p[0] = (char *) calloc(nr * nc, size);
    if (! p[0] ) {
	alloc_error_print(nr * nc, size);
	free(p);
	return NULL;
    }
    for (i=1; i < nr; i++) p[i] = p[i-1] + nc * size;
    return ((void **)p);
}

void free2D(void **p)
{
    free(p[0]);
    free(p);
}

void ***calloc3D(int nr, int nc, int nd, int size)
{
    char ***p, **q, *r;
    int i, j;
    
    p = (char ***)calloc(nr, sizeof(char **));
    if (!p) {
	alloc_error_print(nr, sizeof(char **));
	return NULL;
    }
    q = (char **)calloc(nr * nc,  sizeof(char *));
    if (! q) {
	alloc_error_print(nr * nc,  sizeof(char *));
	free(p);
	return NULL;
    }
    r = (char *)calloc(nr * nc * nd, size);
    if ( ! r){
	alloc_error_print(nr * nc * nd, size);
	free(q);
	free(p);
	return NULL;
    }
    q[0] = r;
    for (i=1;i<nr*nc;i++) {
	q[i] = q[i-1] + nd * size;
    }
    p[0] = q;
    for (i=1;i<nr;i++) {
	p[i] = p[i-1] + nc;
    }
    return ((void ***)p);
}

void free3D(void ***p)
{
    free(**p);
    free(*p);
    free(p);
}

unsigned short **
malloc_matrix_us (long nr, long nc)
{
    unsigned short **m = NULL;
    long i;

    /* Allocate pointers to rows */
    m = (unsigned short **) malloc((size_t) (nr * sizeof(unsigned short *)));
    if (!m) {
        fprintf (stderr, "malloc error in malloc_matrix_us for %ld row pointers.\n", nr);
        m = NULL;
        return m;
    }
    /* Allocate rows and set the pointers to them */
    m[0] = (unsigned short *) malloc((size_t) (nr * nc * sizeof(unsigned short)));
    if (!m[0]) {
        fprintf (stderr, "malloc error in malloc_matrix_us for %ld row with %ld columns.\n", nr, nc);
        m[0] = NULL;
        free (m);
        m = NULL;
        return m;
    }
    for (i = 1; i < nr; i++) m[i] = m[i-1] + nc;

    return m;
}

void
free_matrix_us (unsigned short **m) 
{
    free (m[0]);
    free (m);
    return;
}

