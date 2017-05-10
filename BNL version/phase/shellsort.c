/* Sorting functions */
int
numcmp (char *s1, char *s2)
{
    unsigned short us1, us2;

    us1 = atoi(s1);
    us2 = atoi(s2);
    if (us1 < us2) {
        return -1;
    } else if (us1 > us2) {
        return 1;
    } else {
        return 0;
    }
}

void
shellsort (unsigned long n, unsigned int a[])
{
    unsigned long i, j, v, inc = 1;

    do {
        inc *=3;
        inc++;
    } while (inc <= n);

    do {
        inc /= 3;
        for (i = inc; i < n; i++) {
            v = a[i];
            j = i;
            while (a[j-inc] > v) {
                a[j] = a[j-inc];
                j -= inc;
                if (j < inc) break;
            }
            a[j] = v;
        }
    } while (inc > 1);
    return;
}
