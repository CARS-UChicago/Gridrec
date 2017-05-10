#include <stdio.h>
#include <sys/time.h>

double dt(char mode)
{
  static double t_init;
  double t_now, delta_t;
  struct timeval tv;

  gettimeofday(&tv,0); 
  t_now = (double)tv.tv_sec+(double)tv.tv_usec*1e-6;
  delta_t=0.0;

  switch (mode)  
    {
    case 'i': t_init = t_now;
        break;
    case 'e': delta_t = t_now - t_init; 
        t_init = t_now;
        break;
    case 'r': delta_t = t_now - t_init;
        break;
    default: printf("Invalid timing mode\n");
    }

  return delta_t;
}
