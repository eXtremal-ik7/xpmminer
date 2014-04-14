#include "system.h"
#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include <stdlib.h>
#include <string.h>
#include "config.h"

timeMark getTimeMark()
{
  timeMark mark;  
#ifdef WIN32
  LARGE_INTEGER win32Mark;
  QueryPerformanceCounter(&win32Mark);
  mark.mark = win32Mark.QuadPart;
#else
  struct timeval t;
  gettimeofday(&t, 0);
  mark.mark = t.tv_sec * 1000000 + t.tv_usec;
#endif
  return mark;  
}


uint64_t usDiff(timeMark first, timeMark second)
{
#ifdef WIN32
  LARGE_INTEGER win32Frequency;
  QueryPerformanceFrequency(&win32Frequency);
  return (uint64_t)((second.mark - first.mark) /
                    (double)win32Frequency.QuadPart * 1000000.0);
  return 0;
#else
  uint64_t realSecond = (second.mark >= first.mark) ? second.mark :
    second.mark + 24*3600*(uint64_t)1000000;
  return realSecond - first.mark;    
#endif
}

const char *installPrefix()
{
  return CMAKE_INSTALL_PREFIX;
}

const char *buildPath(PathTy type, const char *fileName)
{
#ifdef WIN32
  char Prefix[MAX_PATH];
  DWORD pathSize =
  GetModuleFileName(GetModuleHandle(NULL), Prefix, sizeof(Prefix));
  
  char *p = Prefix + pathSize;
  while (p > Prefix) {
    if (*p == '\\') {
      *p = 0;
      break;
    }
    p--;
  }
  
  const char *dir = "\\";
#else
  // TODO: сделать директорию установки настраиваемой
  const char Prefix[] = CMAKE_INSTALL_PREFIX;
  
  const char *dir;
  switch (type) {
    case PtExecutable :
      dir = "/bin/";
      break;
    case PtLibrary :
      dir = "/lib/";
      break;
    case PtData :
      dir = "/share/xpmminer/";
      break;
  }
#endif  
  char *path = (char*)
  malloc(sizeof(Prefix)-1 + strlen(dir) + strlen(fileName) + 1);
  strcpy(path, Prefix);
  strcat(path, dir);
  strcat(path, fileName);
  return (const char*)path;
  
}