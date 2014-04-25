#ifndef __DEBUG_H_
#define __DEBUG_H_

#include "stdio.h"
#include "stdarg.h"
#include <time.h>
#include <ncurses.h>

static void dbgPrint(const char *fmt, ...)
{
#ifndef NDEBUG
  va_list arguments;
  va_start(arguments, fmt);
  vfprintf(stderr, fmt, arguments);
#endif
}

static void logFormattedWrite(void *log, const char *fmt, ...)
{
  char buffer[80];       
  time_t rawtime;
  tm *timeinfo;          

  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,80,"%d-%m-%Y %H:%M:%S",timeinfo);   
        
  va_list arguments;
  va_start(arguments, fmt);
  wprintw((WINDOW*)log, "%s: ", buffer);
  vwprintw((WINDOW*)log, fmt, arguments);
  wprintw((WINDOW*)log, "\n");
}

#endif //__DEBUG_H_
