#include "system.h"
#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#include <stdlib.h>
#include <string.h>
#include "config.h"

timeMark getTimeMark()
{
  return std::chrono::steady_clock::now();
}


uint64_t usDiff(timeMark first, timeMark second)
{
  return std::chrono::duration_cast<std::chrono::microseconds>(second-first).count();
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

void xsleep(unsigned seconds)
{
#ifndef WIN32
    sleep(seconds);
#else
    Sleep(seconds*1000);
#endif   
}
