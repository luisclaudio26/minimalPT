#ifndef LOG_H
#define LOG_H

#include <cstdio>

#define ERROR(s) { printf("\nError thrown at %s, line %d: %s\n", __FILE__, __LINE__, s); }

#endif
