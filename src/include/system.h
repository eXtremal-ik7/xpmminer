#include <stdint.h>
#include <chrono>
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (__cplusplus < 201103L)
#define steady_clock monotonic_clock
#endif  

enum PathTy {
  PtExecutable = 0,
  PtLibrary,
  PtData
};

typedef std::chrono::time_point<std::chrono::steady_clock> timeMark;

timeMark getTimeMark();
uint64_t usDiff(timeMark first, timeMark second);

const char *installPrefix();
const char *buildPath(PathTy type, const char *fileName);
void xsleep(unsigned seconds);
