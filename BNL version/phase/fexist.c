#include <sys/stat.h>
int
fexist (char *filename) 
{
    /* Check for existence of a file */
    struct stat stbuf;
    return (stat (filename, &stbuf) == 0);
}
