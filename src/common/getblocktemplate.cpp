#include "getblocktemplate.h"
#include "primecoin.h"

#include <unistd.h>
#include <string>

class MutexLocker {
private:
  pthread_mutex_t *_mutex;
  
public:
  MutexLocker(pthread_mutex_t *mutex) : _mutex(mutex) { pthread_mutex_lock(_mutex); }
  ~MutexLocker() { pthread_mutex_unlock(_mutex); }
};

struct _cbscript_t {
  char *data;
  size_t sz;
};

size_t GetBlockTemplateContext::curlWriteCallback(void *ptr,
                                                  size_t size,
                                                  size_t nmemb,
                                                  GetBlockTemplateContext *ctx)
{
  std::string S((const char*)ptr, size*nmemb);
  ctx->updateWork(S);
  return size*nmemb;
}

void *GetBlockTemplateContext::queryWorkThreadProc(void *arg)
{
  ((GetBlockTemplateContext*)arg)->queryWork();
}

const char *set_b58addr(const char *arg, _cbscript_t *p)
{
  size_t scriptsz = blkmk_address_to_script(NULL, 0, arg);
  if (!scriptsz)
    return "Invalid address";
  
  char *script = (char*)malloc(scriptsz);
  if (blkmk_address_to_script(script, scriptsz, arg) != scriptsz) {
    free(script);
    return "Failed to convert address to script";
  }
  
  p->data = script;
  p->sz = scriptsz;
  return 0;
}

void GetBlockTemplateContext::updateWork(const std::string &S)
{
  MutexLocker locker(&_mutex);
  
  _blockTemplateExists = false;
  const char *error;  
  json_error_t jsonError;
  json_t *response = json_loads(S.c_str(), 0, &jsonError);
  if (!response) {
    fprintf(stderr, " * Error: getblocktemplate response JSON parsing error\n");
    return;
  }
  
  if (_blockTemplate)
    blktmpl_free(_blockTemplate);
  for (unsigned i = 0; i < _blocksNum; i++) {
    if (_blocks[i]) {
      blktmpl_free(_blocks[i]);
      _blocks[i] = 0;
    }
  }
  
  _blockTemplate = blktmpl_create();
  
  error = blktmpl_add_jansson(_blockTemplate, response, time(0));
  json_delete(response);
  if (error) {
    fprintf(stderr, " * Error: %s\n", error);
    return;
  }

  for (unsigned i = 0; i < _blocksNum; i++) {
    unsigned char data[80];
    _cbscript_t opt_coinbase_script;    
     bool newcb;    
   
    error = set_b58addr(_wallet, &opt_coinbase_script);
    if (error) {
      fprintf(stderr, " * Error: %s\n", error);
      return;
    }
    
    _blocks[i] = blktmpl_duplicate(_blockTemplate);
    uint64_t num = 
      blkmk_init_generation2(_blocks[i],
                             opt_coinbase_script.data,
                             opt_coinbase_script.sz,
                             &newcb);
      
    uint32_t extraNonce = (_extraNonce << 16) + i;
    ssize_t coinbaseAppendResult = blkmk_append_coinbase_safe(_blocks[i], &extraNonce, sizeof(extraNonce));
    if (coinbaseAppendResult < 0) {
      fprintf(stderr, " * Error: cannot add extra nonce (error %i)\n", (int)coinbaseAppendResult);
      return; 
    }
      
    size_t dataSize = blkmk_get_data(_blocks[i], data, sizeof(data), time(0), NULL, &_dataId);
    if (!(dataSize >= 76 && dataSize <= sizeof(data))) {
      fprintf(stderr, " * Error: getblocktemplate response decoding error (invalid size)\n");
      return;
    }      
  }
 
  _blockTemplateExists = true;
}

void GetBlockTemplateContext::queryWork()
{
  // JSON-request creation
  char *request;
  {
    blktemplate_t *t = blktmpl_create();
    json_t *jsonRequest = blktmpl_request_jansson(blktmpl_addcaps(t), 0);
    request = json_dumps(jsonRequest, JSON_INDENT(2));
    json_delete(jsonRequest);  
    blktmpl_free(t);
  }
  
  curl_slist *header = curl_slist_append(0, "User-Agent: xpmminer");
  curl_slist_append(header, "Content-Type: application/json");
  
  CURL *curl = curl_easy_init();
  curl_easy_setopt(curl, CURLOPT_URL, _url);
  curl_easy_setopt(curl, CURLOPT_POST, 1L);
  curl_easy_setopt(curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
  curl_easy_setopt(curl, CURLOPT_USERNAME, _user);
  curl_easy_setopt(curl, CURLOPT_PASSWORD, _password);
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriteCallback);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);  
  curl_easy_setopt(curl, CURLOPT_POSTFIELDS, request);
  
  while (1) {
    if (curl_easy_perform(curl) != CURLE_OK)
      fprintf(stderr, "curl_easy_perform ERROR\n");
    sleep(_timeout);
  }
}

GetBlockTemplateContext::GetBlockTemplateContext(const char *url,
                                                 const char *user,
                                                 const char *password,
                                                 const char *wallet,
                                                 unsigned timeout,
                                                 unsigned blocksNum,
                                                 unsigned extraNonce) :
  _url(url), _user(user), _password(password), _wallet(wallet),
  _timeout(timeout), _blocksNum(blocksNum), _extraNonce(extraNonce),
  _blockTemplateExists(false), _blockTemplate(0)
{
  pthread_mutex_init(&_mutex, 0);
  _blocks = (blktemplate_t**)calloc(sizeof(blktemplate_t*), blocksNum);
  printf("workerId=%u\n", extraNonce);
}

void GetBlockTemplateContext::run()
{
  pthread_t thread;
  pthread_create(&thread, 0, queryWorkThreadProc, this);
}

blktemplate_t *GetBlockTemplateContext::get(unsigned blockIdx, blktemplate_t *old, unsigned *dataId, bool *hasChanged)
{
  MutexLocker locker(&_mutex);
  if (old && old->height == _blocks[blockIdx]->height) {
    *hasChanged = false;
    return old;
  } else {
    *hasChanged = true;
    if (old)
      blktmpl_free(old);
    *dataId = _dataId;
    return _blockTemplateExists ? blktmpl_duplicate(_blocks[blockIdx]) : 0;
  }
}


size_t SubmitContext::curlWriteCallback(void *ptr,
                                        size_t size,
                                        size_t nmemb,
                                        void *arg)
{
  std::string S((const char*)ptr, size*nmemb);
  fprintf(stderr, "submit response: %s\n", S.c_str());
  return size*nmemb;
}


SubmitContext::SubmitContext()
{
  curl_slist *header = curl_slist_append(0, "User-Agent: xpmminer");
  curl_slist_append(header, "Content-Type: application/json");
  
  curl = curl_easy_init();
  curl_easy_setopt(curl, CURLOPT_URL, "http://127.0.0.1:9999");
  curl_easy_setopt(curl, CURLOPT_POST, 1L);
  curl_easy_setopt(curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
  curl_easy_setopt(curl, CURLOPT_USERNAME, "user");
  curl_easy_setopt(curl, CURLOPT_PASSWORD, "12345678");
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriteCallback);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);  
}

void SubmitContext::submitBlock(blktemplate_t *blockTemplate,
                                const PrimecoinBlockHeader &header,
                                unsigned dataId)
{
  json_t *jsonBlock =
    blkmk_submit_jansson(blockTemplate,
                         (unsigned char*)&header,
                         dataId,
                         __builtin_bswap32(header.nonce),
                         header.multiplier,
                         header.multiplier[0]+1);
  char *request = json_dumps(jsonBlock, JSON_INDENT(2));
  json_delete(jsonBlock);
  
  fprintf(stderr, "submit request: %s\n", request);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDS, request);
  if (curl_easy_perform(curl) != CURLE_OK)
    fprintf(stderr, "curl_easy_perform ERROR\n");
  
  free(request);
}
