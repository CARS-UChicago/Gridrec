#include <stdlib.h>

main()
{
  long n,m;
  char *p;

  n = 0;
  for(;;) {
    m = 1;
    while ((p = malloc(n+m)) != NULL) {
      m = m + m;
      free(p);
    }
    m /= 2;
    n += m;
    printf ("%ld %lx\n",n,n);
    if (m <= 1) break;
  }
  printf("%ld %lx\n", n,n);
}
