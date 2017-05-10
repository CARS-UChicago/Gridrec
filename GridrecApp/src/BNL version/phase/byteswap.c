#include <stdio.h>
#include <stdlib.h>


long
byteswap (long x, int nbytes) {
    long y = 0;
    unsigned long maskl, maskr;
    int i, n, nd2;

    if (nbytes != (nbytes/2)*2) {
        fprintf (stderr, "Fatal error:"
                " function byteswap only works for an even number of bytes.\n");
        exit(2);
    }

    n = nbytes;
    nd2 = nbytes / 2;
    maskr = ~(~0 << 8);
    maskl = maskr << 8*(nbytes-1);
    
    for (i = 0; i < nd2; i++, n-=2) {
        int yr = 0, yl = 0;
        yr = (x & maskr) << 8*(n-1);
        yl = (x & maskl) >> 8*(n-1);
        y |= yr;
        y |= yl;
        maskr <<= 8;
        maskl >>= 8;
    }
    return y;
}
