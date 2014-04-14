#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

enum PathTy {
  PtExecutable = 0,
  PtLibrary,
  PtData
};
  
struct timeMark {
  uint64_t mark;
};

typedef struct timeMark timeMark;

timeMark getTimeMark();
uint64_t usDiff(timeMark first, timeMark second);

const char *installPrefix();
const char *buildPath(PathTy type, const char *fileName);

#ifdef __cplusplus
}
#endif
