#include "CSieveOfEratosthenesL1Ext.h"
#include "getblocktemplate.h"
#include "primecoin.h"
#include "system.h"
#include "rippedFromHp.h"
#include "Debug.h"

#include <ncurses.h>

#include <getopt.h>
#include <memory>
#include <set>
#include <stdlib.h>

unsigned gDebug = 0;
int gExtensionsNum = 9;
int gPrimorial = 19;
int gSieveSize = CSieveOfEratosthenesL1Ext::L1CacheSize * 10;
int gWeaveDepth = 8192;
int gThreadsNum = 1;
int extraNonce = 0;

static const char *gWallet = 0;
static const char *gUrl = "127.0.0.1:9912";
static const char *gUserName = 0;
static const char *gPassword = 0;


enum CmdLineOptions {
  clDebug = 0,
  clThreadsNum,
  clBenchmark,
  clExtensionsNum,
  clPrimorial,
  clSieveSize,
  clWeaveDepth,
  clUrl,
  clUser,
  clPass,
  clWallet,
  clWorkerId,
  clHelp,
  clOptionLast,
  clOptionsNum
};

void initCmdLineOptions(option *options)
{
  options[clDebug] = {"debug", no_argument, 0, 0};
  options[clThreadsNum] = {"threads", required_argument, &gThreadsNum, 0};
  options[clBenchmark] = {"benchmark", no_argument, 0, 'b'};
  options[clExtensionsNum] = {"extensions-num", required_argument, &gExtensionsNum, 0};  
  options[clPrimorial] = {"primorial", required_argument, &gPrimorial, 0};
  options[clSieveSize] = {"sieve-size", required_argument, &gSieveSize, 0};
  options[clWeaveDepth] = {"weave-depth", required_argument, &gWeaveDepth, 0};
  options[clUrl] = {"url", required_argument, 0, 'o'};
  options[clUser] = {"user", required_argument, 0, 'u'};
  options[clPass] = {"pass", required_argument, 0, 'p'};
  options[clWallet] = {"wallet", required_argument, 0, 'w'};
  options[clWorkerId] = {"worker-id", required_argument, &extraNonce, 0};  
  options[clHelp] = {"help", no_argument, 0, 'h'};
  options[clOptionLast] = {0, 0, 0, 0};
}

void sieveL1ExtBenchmark(double difficulty)
{
  const unsigned HpSieveExtensionsNum = 9;
  const unsigned HpSieveL1CacheSize = 220400;        
        
  PrimeSource primeSource(1000000, gWeaveDepth);
  CPrimalityTestParams testParams(bitsFromDifficulty(difficulty));
  
  PrimecoinBlockHeader header;
  mpz_class blockHeaderHash;
  mpz_class primorial;  
  mpz_class fixedMultiplier;
  generateRandomHeader(&header, difficulty);
  if (!updateBlock(&header, blockHeaderHash, primeSource, testParams)) {
    fprintf(stderr, "sieveBenchmark error: cannot find block with correct SHA256 hash\n");
    exit(1);
  }
  
  PrimorialFast(gPrimorial, primorial, primeSource);
  fixedMultiplier = primorial * blockHeaderHash;
  for (unsigned i = 1; i <= 8; i++) {
    unsigned sieveSize = CSieveOfEratosthenesL1Ext::L1CacheSize*2*i;
    unsigned realSieveSize =
      sieveSize + sieveSize/2*CSieveOfEratosthenesL1Ext::ExtensionsNum;
      
      // Fast checking CSieveOfEratosthenesL1Ext output
    {
      std::unique_ptr<CSieveOfEratosthenesL1Ext> l1ext(
        new CSieveOfEratosthenesL1Ext(primeSource));
      std::unique_ptr<CSieveOfEratosthenesHp> hp(
        new CSieveOfEratosthenesHp(primeSource));
      l1ext->reset(sieveSize, (unsigned)difficulty, gWeaveDepth, fixedMultiplier);      
      hp->Reset(sieveSize, gWeaveDepth, HpSieveExtensionsNum, HpSieveL1CacheSize,
                header.bits, blockHeaderHash, primorial, 0);
      l1ext->Weave();
      hp->Weave();
      if (l1ext->GetCandidateCount() != hp->GetCandidateCount() || !l1ext->fastSelfTest()) {
        fprintf(stderr, "CSieveOfEratosthenesL1Ext weave failed!\n");
        exit(1);
      }
    }

    unsigned limit;
    if (sieveSize < 100000) limit = 60;
    else if (sieveSize < 500000) limit = 30;
    else limit = 15;
    
    timeMark beginPoint = getTimeMark();    
    for (unsigned j = 0; j < limit; j++) {
      std::unique_ptr<CSieveOfEratosthenesL1Ext> sieve(
        new CSieveOfEratosthenesL1Ext(primeSource));
      sieve->reset(sieveSize, (unsigned)difficulty, gWeaveDepth, fixedMultiplier);
      sieve->Weave();      
    }
    
    timeMark middlePoint = getTimeMark();        
    for (unsigned j = 0; j < limit; j++) {
      std::unique_ptr<CSieveOfEratosthenesHp> sieve(
        new CSieveOfEratosthenesHp(primeSource));
      sieve->Reset(sieveSize, gWeaveDepth,
                   HpSieveExtensionsNum, HpSieveL1CacheSize,
                   header.bits, blockHeaderHash, primorial, 0);
      sieve->Weave();      
    }
    
    timeMark endPoint = getTimeMark();        
    double tm = usDiff(beginPoint, middlePoint) / (1000.0 * limit);
    double hpTm = usDiff(middlePoint, endPoint) / (1000.0 * limit);
    fprintf(stderr, "(%u, real %u) L1Ext(time %.3lfms speed %.3lf) HP(time %.3lfms speed %.3lf) %.3lf%% faster\n",
            sieveSize, realSieveSize,
            tm, realSieveSize/tm,
            hpTm, realSieveSize/hpTm,
            hpTm/tm*100.0 - 100.0);
  }
}

void fermatTestBenchmark(double difficulty)
{
  PrimeSource primeSource(1000000, gWeaveDepth);
  CPrimalityTestParams testParams(bitsFromDifficulty(difficulty));
  
  PrimecoinBlockHeader header;
  mpz_class blockHeaderHash;
  mpz_class primorial;  
  mpz_class fixedMultiplier;
  generateRandomHeader(&header, difficulty);

  while (updateBlock(&header, blockHeaderHash, primeSource, testParams)) {
    PrimorialFast(gPrimorial, primorial, primeSource);
    fixedMultiplier = primorial * blockHeaderHash;

    std::vector<mpz_class> primesForTest;    
    std::unique_ptr<CSieveOfEratosthenesL1Ext> sieve(
      new CSieveOfEratosthenesL1Ext(primeSource));
    sieve->reset(gSieveSize, (unsigned)difficulty, gWeaveDepth, fixedMultiplier);
    sieve->Weave();
    sieve->resetCandidateIterator();
      
    mpz_class prime;
    unsigned testsNum = 0;
    unsigned multiplier;
    unsigned candidateType;

    while (sieve->GetNextCandidateMultiplier(multiplier, candidateType)) {
      prime = fixedMultiplier;
      prime *= multiplier;
      if (candidateType == PRIME_CHAIN_CUNNINGHAM1 ||
          candidateType == PRIME_CHAIN_BI_TWIN) {
        primesForTest.push_back(prime-1);
      }
      
      if (candidateType == PRIME_CHAIN_CUNNINGHAM2 ||
          candidateType == PRIME_CHAIN_BI_TWIN) {
        primesForTest.push_back(prime+1);
      }      
    }
    
    timeMark beginPoint = getTimeMark();
    for (size_t i = 0; i < primesForTest.size(); i++) {
      unsigned chainLength;
      FermatProbablePrimalityTestFast(primesForTest[i], chainLength, testParams);
    }
    uint64_t testTime = usDiff(beginPoint, getTimeMark());
    
    fprintf(stderr,
            " %u Fermat tests at %.3lfms, speed: %.3lf tests/second\n",
            (unsigned)primesForTest.size(),
            testTime / 1000.0,
            1000000.0 / testTime * primesForTest.size());
  }
}

void benchmark(double difficulty)
{
  unsigned realSieveSize =
    gSieveSize + gSieveSize/2*CSieveOfEratosthenesL1Ext::ExtensionsNum;  
  PrimeSource primeSource(1000000, gWeaveDepth);
  CPrimalityTestParams testParams(bitsFromDifficulty(difficulty));
  
  std::unique_ptr<CSieveOfEratosthenesL1Ext> sieve(
    new CSieveOfEratosthenesL1Ext(primeSource));
  uint64_t foundChains[20];
  memset(foundChains, 0, sizeof(foundChains));
  
  while (true) {
    // Create random primecoin header
    PrimecoinBlockHeader header;
    generateRandomHeader(&header, difficulty);
    
    double totalRoundTime = 0.0;
    double totalSpeed = 0.0;
    unsigned roundsNum = 0;
    mpz_class blockHeaderHash;
    mpz_class primorial;
    while (true) {
      timeMark roundBegin = getTimeMark();
      
      // Create suitable block by nonce searching
      if (!updateBlock(&header, blockHeaderHash, primeSource, testParams))
        break;

      PrimorialFast(gPrimorial, primorial, primeSource);
      
      // Cunningham chain search
      unsigned probableChainLength;
      unsigned testsNum;
      unsigned primeHitsNum;
      MineProbablePrimeChainFast(header,
                                 sieve.get(),
                                 blockHeaderHash,
                                 primorial,
                                 probableChainLength,
                                 testsNum,
                                 primeHitsNum,
                                 testParams,
                                 primeSource,
                                 foundChains);
      
      double roundTime = usDiff(roundBegin, getTimeMark());
      if (gDebug) {
        fprintf(stderr,
                "  * total round time: %.3lf milliseconds, speed: %.3lf M\n",
                roundTime / 1000.0,
                realSieveSize / roundTime);
      }
      
      totalRoundTime += (roundTime / 1000.0);
      totalSpeed += (realSieveSize / roundTime);
      roundsNum++;
      if (!(roundsNum % 16))
        fprintf(stderr, "average speed: %.3lf M\n", totalSpeed / roundsNum);
    }
  }
}


struct MineContext {
  PrimeSource *primeSource;
  GetBlockTemplateContext *gbp;
  SubmitContext *submit;
  unsigned threadIdx;
  uint64_t totalRoundsNum; 
  uint64_t foundChains[20];
  double speed;
  WINDOW *log;
};

void *mine(void *arg)
{
  MineContext *ctx = (MineContext*)arg;
  unsigned realSieveSize =
    gSieveSize + gSieveSize/2*CSieveOfEratosthenesL1Ext::ExtensionsNum;  

  blktemplate_t *workTemplate = 0;
  PrimecoinBlockHeader work;
  unsigned dataId;
  mpz_class blockHeaderHash;
  mpz_class primorial;

  const unsigned checkInterval = 8;
  double roundSizeInGb = checkInterval*realSieveSize / 1000000000.0;
  unsigned roundsNum = 0;    
  
  PrimorialFast(gPrimorial, primorial, *ctx->primeSource);
  std::unique_ptr<CSieveOfEratosthenesL1Ext> sieve(
    new CSieveOfEratosthenesL1Ext(*ctx->primeSource));    
  CPrimalityTestParams testParams(0);  

  timeMark localWorkBegin = getTimeMark();
  while (1) {
    bool hasChanged;
    while ( !(workTemplate = ctx->gbp->get(ctx->threadIdx, workTemplate, &dataId, &hasChanged)) )
      usleep(100);
    
    timeMark roundBegin = getTimeMark();    

    if (hasChanged) {
      work.version = workTemplate->version;
      memcpy(work.hashPrevBlock, workTemplate->prevblk, 32);
      memcpy(work.hashMerkleRoot, workTemplate->_mrklroot, 32);
      work.time = workTemplate->curtime;
      work.bits = *(uint32_t*)workTemplate->diffbits;
      work.nonce = 0;
      testParams.bits = work.bits;
    }

    if (!updateBlock(&work, blockHeaderHash, *ctx->primeSource, testParams))
      continue;
    
    unsigned probableChainLength;
    unsigned testsNum;
    unsigned primeHitsNum;
    if (MineProbablePrimeChainFast(work,
                                   sieve.get(),
                                   blockHeaderHash,
                                   primorial,
                                   probableChainLength,
                                   testsNum,
                                   primeHitsNum,
                                   testParams,
                                   *ctx->primeSource,
                                   ctx->foundChains)) {
      logFormattedWrite(ctx->log, "block found!");
      ctx->submit->submitBlock(workTemplate, work, dataId);
    }
      
    roundsNum++;
    if (roundsNum == checkInterval) {
      timeMark point = getTimeMark();
      uint64_t timeElapsed = usDiff(localWorkBegin, point);
      
      ctx->totalRoundsNum += checkInterval;
      ctx->speed = roundSizeInGb / (timeElapsed / 1000000.0);

      roundsNum = 0;
      localWorkBegin = getTimeMark();
    }
  }
}

void printHelpMessage() 
{
  printf("Opensource primecoin CPU miner, usage:\n");
  printf("  xpmclminer <arguments>\n\n");
  printf("  -h or --help: show this help message\n");
  printf("  -b or --benchmark: run benchmark and exit\n");
  printf("  -o or --url <HostAddress:port>: address of primecoin RPC client, default: %s\n", gUrl);
  printf("  -u or --user <UserName>: user name for primecoin RPC client\n");
  printf("  -p or --pass <Password>: password for primecoin RPC client\n");
  printf("  -w or --wallet: wallet address for coin receiving\n");
  printf("  --debug: show additional mining information\n");
  printf("  --extensions-num <number>: Eratosthenes sieve extensions number (default: %u)\n", gExtensionsNum);
  printf("  --primorial <number>: primorial number (default: %u)\n", gPrimorial);
  printf("  --sieve-size <number>: Eratosthenes sieve size (default: %u)\n", gSieveSize);
  printf("  --weave-depth <number>: Eratosthenes sieve weave depth (default: %u)\n", gWeaveDepth);
  printf("  --worker-id: unique identifier of your worker, used in block creation. ");
  printf("All your rigs must have different worker IDs! (default: current time value)\n");
}

int main(int argc, char **argv)
{
  srand(time(0));  
  blkmk_sha256_impl = sha256;  
  option gOptions[clOptionsNum];
  PrimeSource primeSource(10000000, gWeaveDepth+256);
  
  bool isBenchmark = false;
  int index = 0, c;
  initCmdLineOptions(gOptions);
  const char *platform = "Advanced Micro Devices, Inc.";
  while ((c = getopt_long(argc, argv, "bo:u:p:w:h", gOptions, &index)) != -1) {
    switch (c) {
      case 0 :
        switch (index) {
          case clDebug :
            gDebug = 1;
            break;
          case clExtensionsNum :
            gExtensionsNum = atoi(optarg);
            break;
          case clPrimorial :
            gPrimorial = atoi(optarg);
            break;
          case clSieveSize :
            gSieveSize = atoi(optarg);
            if (gSieveSize > CSieveOfEratosthenesL1Ext::MaxSieveSize) {
              fprintf(stderr, "Sieve size limited by %u, you try launch with %u\n",
                      CSieveOfEratosthenesL1Ext::MaxSieveSize, gSieveSize);
              exit(1);
            }
            break;
          case clWeaveDepth :
            gWeaveDepth = atoi(optarg);
            if (gWeaveDepth > CSieveOfEratosthenesL1Ext::MaxDepth) {
              fprintf(stderr, "Weave depth limited by %u, you try launch with %u\n",
                      CSieveOfEratosthenesL1Ext::MaxDepth, gWeaveDepth);
              exit(1);
            }
            break;
          case clThreadsNum :
            gThreadsNum = atoi(optarg);
            break;
          case clWorkerId :
            extraNonce = atoi(optarg);
            break;            
        }
        break;
          case 'b' :
            isBenchmark = true;
            break;
          case 'o' :
            gUrl = optarg;
            break;
          case 'u' :
            gUserName = optarg;
            break;
          case 'p' :
            gPassword = optarg;
            break;
          case 'w' :
            gWallet = optarg;
            break;
          case 'h' :
            printHelpMessage();
            exit(0);
          case ':' :
            fprintf(stderr, "Error: option %s missing argument\n",
                    gOptions[index].name);
            break;
          case '?' :
            fprintf(stderr, "Error: invalid option %s\n", argv[optind-1]);
            break;
          default :
            break;
    }
  }
  
  if ((!gUserName || !gPassword) && !isBenchmark) {
    fprintf(stderr, "Error: you must specify user name and password\n");
    exit(1);
  }
  
  if (!gWallet && !isBenchmark) {
    fprintf(stderr, "Error: you must specify wallet\n");
    exit(1);
  }
  
  if (extraNonce == 0)
    extraNonce = time(0);

  if (isBenchmark) {
    fermatTestBenchmark(10.5);
    sieveL1ExtBenchmark(10.5);
    benchmark(10.5);
    return 0;
  }
 
  WINDOW *display = initscr();
  WINDOW *log = newwin(30, 160, 3 + gThreadsNum + 12, 0);
  scrollok(log, TRUE); 
 
  GetBlockTemplateContext ctx(log, gUrl, gUserName, gPassword, gWallet, 4, gThreadsNum, extraNonce);
  ctx.run();
  
  MineContext *mineCtx = new MineContext[gThreadsNum];
  for (unsigned i = 0; i < gThreadsNum; i++) {
    pthread_t thread;    
    mineCtx[i].primeSource = &primeSource;
    mineCtx[i].gbp = &ctx;
    mineCtx[i].threadIdx = i;
    mineCtx[i].totalRoundsNum = 0;
    memset(mineCtx[i].foundChains, 0, sizeof(mineCtx->foundChains));
    mineCtx[i].submit = new SubmitContext(log, gUrl, gUserName, gPassword);
    mineCtx[i].log = log;
    pthread_create(&thread, 0, mine, &mineCtx[i]);
  }
  
  unsigned realSieveSize =
    gSieveSize + gSieveSize/2*CSieveOfEratosthenesL1Ext::ExtensionsNum;    
  double sieveSizeInGb = realSieveSize / 1000000000.0;
  timeMark workBeginPoint = getTimeMark();

  {
    time_t rawtime;
    struct tm * timeinfo;          
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,80, "%d-%m-%Y %H:%M:%S", timeinfo);  
    wmove(display, 0, 0);
    wprintw(display, " ** xpmcpuminer started %s %s %s worker %i **\n", buffer, gUrl, gWallet, extraNonce);
  }
  
  refresh();
  while (true) {
    xsleep(5);  
    uint64_t foundChains[MaxChainLength];
    double speed = 0.0;
    double averageSpeed = 0.0;
    memset(foundChains, 0, sizeof(foundChains));
    
    wmove(display, 1, 0);
    wprintw(display, " ** block: %u, difficulty: %.3lf", ctx.getBlockHeight(), ctx.getDifficulty());
    timeMark currentPoint = getTimeMark();    
    uint64_t elapsedTime = usDiff(workBeginPoint, currentPoint);
    for (int i = 0; i < gThreadsNum; i++) {
      for (unsigned chIdx = 1; chIdx < MaxChainLength; chIdx++)
        foundChains[chIdx] += mineCtx[i].foundChains[chIdx];
      
      double threadAvgSpeed = (sieveSizeInGb*mineCtx[i].totalRoundsNum) / (elapsedTime / 1000000.0);
      speed += mineCtx[i].speed;
      averageSpeed += threadAvgSpeed;
      
      wmove(display, i+3, 0);
      wprintw(display, "[%u] %.3lfG, average: %.3lfG", i+1, mineCtx[i].speed, threadAvgSpeed);
    }
  
    wmove(display, gThreadsNum+3, 0);
    wprintw(display, " * speed: %.3lfG, average: %.3lfG\n", speed, averageSpeed);
    unsigned chIdx;
    for (chIdx = 1; chIdx < MaxChainLength && foundChains[chIdx]; chIdx++) {
      wmove(display, gThreadsNum+3 + chIdx+1, 0);
      wprintw(display, "   * chains/%u: %llu %.3lf/sec ",
              chIdx, foundChains[chIdx], foundChains[chIdx] / (elapsedTime / 1000000.0));
      if (chIdx >= 7)
        wprintw(display, "%.3lf/hour ", foundChains[chIdx] / (elapsedTime / 1000000.0) * 3600.0);
    }
    
    wrefresh(display);
    wrefresh(log);
  }
  
  return 0;
}
