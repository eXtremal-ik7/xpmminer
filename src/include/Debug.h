#ifndef __DEBUG_H_
#define __DEBUG_H_

#include "stdio.h"
#include "stdarg.h"

static void dbgPrint(const char *fmt, ...)
{
#ifndef NDEBUG
  va_list arguments;
  va_start(arguments, fmt);
  vfprintf(stderr, fmt, arguments);
#endif
}

#endif //__DEBUG_H_
