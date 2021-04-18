/*
 * sha256.cpp
 *
 *  Created on: 11.01.2014
 *      Author: mad
 */


#include <cstring>
#include <fstream>
#include "sha256.h"
 
const SHA_256::uint32 SHA_256::sha256_k[64] = //UL = uint32
            {0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
             0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
             0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
             0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
             0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
             0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
             0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
             0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
             0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
             0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
             0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
             0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
             0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
             0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
             0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
             0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
 
void SHA_256::transform(const unsigned char *message, unsigned int block_nb)
{
    uint32 w[64];
    uint32 wv[8];
    uint32 t1, t2;
    const unsigned char *sub_block;
    int i;
    int j;
    for (i = 0; i < (int) block_nb; i++) {
        sub_block = message + (i << 6);
        for (j = 0; j < 16; j++) {
            SHA2_PACK32(&sub_block[j << 2], &w[j]);
        }
        for (j = 16; j < 64; j++) {
            w[j] =  SHA256_F4(w[j -  2]) + w[j -  7] + SHA256_F3(w[j - 15]) + w[j - 16];
        }
        for (j = 0; j < 8; j++) {
            wv[j] = m_h[j];
        }
        for (j = 0; j < 64; j++) {
            t1 = wv[7] + SHA256_F2(wv[4]) + SHA2_CH(wv[4], wv[5], wv[6])
                + sha256_k[j] + w[j];
            t2 = SHA256_F1(wv[0]) + SHA2_MAJ(wv[0], wv[1], wv[2]);
            wv[7] = wv[6];
            wv[6] = wv[5];
            wv[5] = wv[4];
            wv[4] = wv[3] + t1;
            wv[3] = wv[2];
            wv[2] = wv[1];
            wv[1] = wv[0];
            wv[0] = t1 + t2;
        }
        for (j = 0; j < 8; j++) {
            m_h[j] += wv[j];
        }
    }
}
 
void SHA_256::init()
{
    m_h[0] = 0x6a09e667;
    m_h[1] = 0xbb67ae85;
    m_h[2] = 0x3c6ef372;
    m_h[3] = 0xa54ff53a;
    m_h[4] = 0x510e527f;
    m_h[5] = 0x9b05688c;
    m_h[6] = 0x1f83d9ab;
    m_h[7] = 0x5be0cd19;
    m_len = 0;
    m_tot_len = 0;
}
 
void SHA_256::update(const unsigned char *message, unsigned int len)
{
    unsigned int block_nb;
    unsigned int new_len, rem_len, tmp_len;
    const unsigned char *shifted_message;
    tmp_len = SHA224_256_BLOCK_SIZE - m_len;
    rem_len = len < tmp_len ? len : tmp_len;
    memcpy(&m_block[m_len], message, rem_len);
    if (m_len + len < SHA224_256_BLOCK_SIZE) {
        m_len += len;
        return;
    }
    new_len = len - rem_len;
    block_nb = new_len / SHA224_256_BLOCK_SIZE;
    shifted_message = message + rem_len;
    transform(m_block, 1);
    transform(shifted_message, block_nb);
    rem_len = new_len % SHA224_256_BLOCK_SIZE;
    memcpy(m_block, &shifted_message[block_nb << 6], rem_len);
    m_len = rem_len;
    m_tot_len += (block_nb + 1) << 6;
}
 
void SHA_256::final(unsigned char *digest)
{
    unsigned int block_nb;
    unsigned int pm_len;
    unsigned int len_b;
    int i;
    block_nb = (1 + ((SHA224_256_BLOCK_SIZE - 9)
                     < (m_len % SHA224_256_BLOCK_SIZE)));
    len_b = (m_tot_len + m_len) << 3;
    pm_len = block_nb << 6;
    memset(m_block + m_len, 0, pm_len - m_len);
    m_block[m_len] = 0x80;
    SHA2_UNPACK32(len_b, m_block + pm_len - 4);
    transform(m_block, block_nb);
    for (i = 0 ; i < 8; i++) {
        SHA2_UNPACK32(m_h[i], &digest[i << 2]);
    }
}
 
std::string sha256_str(std::string input)
{
    unsigned char digest[SHA_256::DIGEST_SIZE];
    memset(digest,0,SHA_256::DIGEST_SIZE);
 
    SHA_256 ctx = SHA_256();
    ctx.init();
    ctx.update( (unsigned char*)input.c_str(), input.length());
    ctx.final(digest);
 
    char buf[2*SHA_256::DIGEST_SIZE+1];
    buf[2*SHA_256::DIGEST_SIZE] = 0;
    for (unsigned i = 0; i < SHA_256::DIGEST_SIZE; i++)
        sprintf(buf+i*2, "%02x", digest[i]);
    return std::string(buf);
}

const uint32_t SHA256_K[64] = {
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
  0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
  0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
  0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
  0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
  0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
  0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
  0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
  0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

uint32_t sha2_pack(uint32_t val) {
  
  return ((val & 0xFF) << 24) |
         ((val & 0xFF00) << 8) |
         ((val & 0xFF0000) >> 8) |
         ((val & 0xFF000000) >> 24);
  
}

#define Zrotr(x,y) ((x<<y) | (x>>(sizeof(x)*8-y)))
#define ZR25(n) ((Zrotr((n), 25) ^ Zrotr((n), 14) ^ ((n) >> 3U)))
#define ZR15(n) ((Zrotr((n), 15) ^ Zrotr((n), 13) ^ ((n) >> 10U)))
#define ZR26(n) ((Zrotr((n), 26) ^ Zrotr((n), 21) ^ Zrotr((n), 7)))
#define ZR30(n) ((Zrotr((n), 30) ^ Zrotr((n), 19) ^ Zrotr((n), 10)))
#define Ch(x, y, z) (z ^ (x & (y ^ z)))
#define Ma(x, y, z) ((x & z) | (y & (x | z)))

void sha256Export(const uint32_t *msg, uint32_t *s,
                  uint32_t *w,
                  uint32_t *temp1Data, uint32_t *temp2Data,
                  uint32_t *new1Data, uint32_t *new2Data)
{
  for(int i = 0; i < 16; ++i)
    w[i] = msg[i];

  for(int i = 16; i < 64; ++i){
    const uint32_t s0 = ZR25(w[i-15]);
    const uint32_t s1 = ZR15(w[i-2]);
    w[i] = w[i-16] + s0 + w[i-7] + s1;    
  }
 
  uint32_t a = s[0];
  uint32_t b = s[1];
  uint32_t c = s[2];
  uint32_t d = s[3];
  uint32_t e = s[4];
  uint32_t f = s[5];
  uint32_t g = s[6];
  uint32_t h = s[7];
  
  for(int i = 0; i < 64; ++i){
    const uint32_t S1 = ZR26(e);
    const uint32_t ch = Ch(e, f, g);
     
    const uint32_t temp1 = h + S1 + ch + SHA256_K[i] + w[i];
    const uint32_t S0 = ZR30(a);
    const uint32_t maj = Ma(a, b, c);
    const uint32_t temp2 = S0 + maj;
    
    h = g;
    g = f;
    f = e;
    e = d + temp1;
    d = c;
    c = b;
    b = a;
    a = temp1 + temp2;
    
    temp1Data[i] = temp1;
    temp2Data[i] = temp2;
    new1Data[i] = a;
    new2Data[i] = e;
  }
  
  s[0] += a;
  s[1] += b;
  s[2] += c;
  s[3] += d;
  s[4] += e;
  s[5] += f;
  s[6] += g;
  s[7] += h; 
}

void precalcSHA256(const void *block,
                   uint32_t *midstate,
                   sha256precalcData *data)
{
  SHA_256 sha;
  sha.init();
  sha.transform((const unsigned char*)block, 1u);
  for(int i = 0; i < 8; ++i)
    midstate[i] = sha.m_h[i];
  
  uint32_t *p = ((uint32_t*)block)+16;   
  uint32_t msg_merkle = sha2_pack(p[0]);
  uint32_t msg_time = sha2_pack(p[1]);
  uint32_t msg_bits = sha2_pack(p[2]);
  
  uint32_t msg[16];
  msg[0] = msg_merkle;
  msg[1] = msg_time;
  msg[2] = msg_bits;
  msg[4] = sha2_pack(0x80);
  for(int i = 5; i < 15; ++i)
    msg[i] = 0;
  msg[15] = 640;  
  
  uint32_t out[8];
  memcpy(out, midstate, 8*4);
  uint32_t W[64];       // 2 [16-17]
  uint32_t temp1[64];   // 3
  uint32_t temp2[64];   // 4
  uint32_t new1[64];    // 3
  uint32_t new2[64];    // 3
  sha256Export(msg, out, W, temp1, temp2, new1, new2);
  
  data->merkle = msg_merkle;
  data->time = msg_time;
  data->nbits = msg_bits;
  data->W0 = W[16];
  data->W1 = W[17];
  data->new1_0 = new1[0];
  data->new1_1 = new1[1];
  data->new1_2 = new1[2];
  data->new2_0 = new2[0];
  data->new2_1 = new2[1];
  data->new2_2 = new2[2];
  data->temp2_3 = temp2[3];
}
