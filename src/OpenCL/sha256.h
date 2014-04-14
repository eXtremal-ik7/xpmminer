__constant uint ES[2] = { 0x00FF00FF, 0xFF00FF00 };
// 
// __constant uint K256[64] = {
//   0x428a2f98UL, 0x71374491UL, 0xb5c0fbcfUL, 0xe9b5dba5UL,
//   0x3956c25bUL, 0x59f111f1UL, 0x923f82a4UL, 0xab1c5ed5UL,
//   0xd807aa98UL, 0x12835b01UL, 0x243185beUL, 0x550c7dc3UL,
//   0x72be5d74UL, 0x80deb1feUL, 0x9bdc06a7UL, 0xc19bf174UL,
//   0xe49b69c1UL, 0xefbe4786UL, 0x0fc19dc6UL, 0x240ca1ccUL,
//   0x2de92c6fUL, 0x4a7484aaUL, 0x5cb0a9dcUL, 0x76f988daUL,
//   0x983e5152UL, 0xa831c66dUL, 0xb00327c8UL, 0xbf597fc7UL,
//   0xc6e00bf3UL, 0xd5a79147UL, 0x06ca6351UL, 0x14292967UL,
//   0x27b70a85UL, 0x2e1b2138UL, 0x4d2c6dfcUL, 0x53380d13UL,
//   0x650a7354UL, 0x766a0abbUL, 0x81c2c92eUL, 0x92722c85UL,
//   0xa2bfe8a1UL, 0xa81a664bUL, 0xc24b8b70UL, 0xc76c51a3UL,
//   0xd192e819UL, 0xd6990624UL, 0xf40e3585UL, 0x106aa070UL,
//   0x19a4c116UL, 0x1e376c08UL, 0x2748774cUL, 0x34b0bcb5UL,
//   0x391c0cb3UL, 0x4ed8aa4aUL, 0x5b9cca4fUL, 0x682e6ff3UL,
//   0x748f82eeUL, 0x78a5636fUL, 0x84c87814UL, 0x8cc70208UL,
//   0x90befffaUL, 0xa4506cebUL, 0xbef9a3f7UL, 0xc67178f2UL
// };
// 
// #define EndianSwap(n) (rotate(n & ES[0], 24U)|rotate(n & ES[1], 8U))
// 
// /* Shift-right (used in SHA-256, SHA-384, and SHA-512): */
// #define R(b,x)    ((x) >> (b))
// /* 32-bit Rotate-right (used in SHA-256): */
// #define S32(b,x)  (((x) >> (b)) | ((x) << (32 - (b))))
// 
// /* Two of six logical functions used in SHA-256, SHA-384, and SHA-512: */
// #define Ch(x,y,z) (((x) & (y)) ^ ((~(x)) & (z)))
// #define Maj(x,y,z)  (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
// 
// /* Four of six logical functions used in SHA-256: */
// #define Sigma0_256(x) (S32(2,  (x)) ^ S32(13, (x)) ^ S32(22, (x)))
// #define Sigma1_256(x) (S32(6,  (x)) ^ S32(11, (x)) ^ S32(25, (x)))
// #define sigma0_256(x) (S32(7,  (x)) ^ S32(18, (x)) ^ R(3 ,   (x)))
// #define sigma1_256(x) (S32(17, (x)) ^ S32(19, (x)) ^ R(10,   (x)))
// 
// #define ROUND256_0_TO_15(a,b,c,d,e,f,g,h, data) \
//   T1 = (h) + Sigma1_256(e) + Ch((e), (f), (g)) + \
//        K256[j] + (W256[j] = data); \
//   (d) += T1; \
//   (h) = T1 + Sigma0_256(a) + Maj((a), (b), (c)); \
//   j++
// 
// #define ROUND256(a,b,c,d,e,f,g,h) \
//   s0 = W256[(j+1)&0x0f]; \
//   s0 = sigma0_256(s0); \
//   s1 = W256[(j+14)&0x0f]; \
//   s1 = sigma1_256(s1); \
//   T1 = (h) + Sigma1_256(e) + Ch((e), (f), (g)) + K256[j] + \
//        (W256[j&0x0f] += s1 + W256[(j+9)&0x0f] + s0); \
//   (d) += T1; \
//   (h) = T1 + Sigma0_256(a) + Maj((a), (b), (c)); \
//   j++  
//   
//   
//   
// void sha256(uint4*restrict state0,uint4*restrict state1, const uint4 block0, const uint4 block1, const uint4 block2, const uint4 block3)
// {
//   uint4 S0 = *state0;
//   uint4 S1 = *state1;
//   
// #define a S0.x
// #define b S0.y
// #define c S0.z
// #define d S0.w
// #define e S1.x
// #define f S1.y
// #define g S1.z
// #define h S1.w
//   
//   uint W256[16], T1, s0, s1;
//   
//   uint j = 0;
//   ROUND256_0_TO_15(a,b,c,d,e,f,g,h, block0.x);
//   ROUND256_0_TO_15(h,a,b,c,d,e,f,g, block0.y);
//   ROUND256_0_TO_15(g,h,a,b,c,d,e,f, block0.z);
//   ROUND256_0_TO_15(f,g,h,a,b,c,d,e, block0.w);
//   ROUND256_0_TO_15(e,f,g,h,a,b,c,d, block1.x);
//   ROUND256_0_TO_15(d,e,f,g,h,a,b,c, block1.y);
//   ROUND256_0_TO_15(c,d,e,f,g,h,a,b, block1.z);
//   ROUND256_0_TO_15(b,c,d,e,f,g,h,a, block1.w);
// //   printf("<-- %08X %08X %08X %08X %08X %08X %08X %08X\n", a, b, c, d, e, f, g, h);        
//   ROUND256_0_TO_15(a,b,c,d,e,f,g,h, block2.x);
//   ROUND256_0_TO_15(h,a,b,c,d,e,f,g, block2.y);
//   ROUND256_0_TO_15(g,h,a,b,c,d,e,f, block2.z);
//   ROUND256_0_TO_15(f,g,h,a,b,c,d,e, block2.w);
//   ROUND256_0_TO_15(e,f,g,h,a,b,c,d, block3.x);
//   ROUND256_0_TO_15(d,e,f,g,h,a,b,c, block3.y);
//   ROUND256_0_TO_15(c,d,e,f,g,h,a,b, block3.z);
//   ROUND256_0_TO_15(b,c,d,e,f,g,h,a, block3.w);  
// //   printf("<-- %08X %08X %08X %08X %08X %08X %08X %08X\n", a, b, c, d, e, f, g, h);      
//   
//   do {
//     ROUND256(a,b,c,d,e,f,g,h);
//     ROUND256(h,a,b,c,d,e,f,g);
//     ROUND256(g,h,a,b,c,d,e,f);
//     ROUND256(f,g,h,a,b,c,d,e);
//     ROUND256(e,f,g,h,a,b,c,d);
//     ROUND256(d,e,f,g,h,a,b,c);
//     ROUND256(c,d,e,f,g,h,a,b);
//     ROUND256(b,c,d,e,f,g,h,a);
// //     printf("<!!-- %08X %08X %08X %08X %08X %08X %08X %08X\n", a, b, c, d, e, f, g, h);          
//   } while (j < 64);
// #undef a
// #undef b
// #undef c
// #undef d
// #undef e
// #undef f
// #undef g
// #undef h
//   
//   *state0 += S0;
//   *state1 += S1;  
// }

__constant uint K[] = {
  0x428a2f98U,
  0x71374491U,
  0xb5c0fbcfU,
  0xe9b5dba5U,
  0x3956c25bU,
  0x59f111f1U,
  0x923f82a4U,
  0xab1c5ed5U,
  0xd807aa98U,
  0x12835b01U,
  0x243185beU, // 10
  0x550c7dc3U,
  0x72be5d74U,
  0x80deb1feU,
  0x9bdc06a7U,
  0xe49b69c1U,
  0xefbe4786U,
  0x0fc19dc6U,
  0x240ca1ccU,
  0x2de92c6fU,
  0x4a7484aaU, // 20
  0x5cb0a9dcU,
  0x76f988daU,
  0x983e5152U,
  0xa831c66dU,
  0xb00327c8U,
  0xbf597fc7U,
  0xc6e00bf3U,
  0xd5a79147U,
  0x06ca6351U,
  0x14292967U, // 30
  0x27b70a85U,
  0x2e1b2138U,
  0x4d2c6dfcU,
  0x53380d13U,
  0x650a7354U,
  0x766a0abbU,
  0x81c2c92eU,
  0x92722c85U,
  0xa2bfe8a1U,
  0xa81a664bU, // 40
  0xc24b8b70U,
  0xc76c51a3U,
  0xd192e819U,
  0xd6990624U,
  0xf40e3585U,
  0x106aa070U,
  0x19a4c116U,
  0x1e376c08U,
  0x2748774cU,
  0x34b0bcb5U, // 50
  0x391c0cb3U,
  0x4ed8aa4aU,
  0x5b9cca4fU,
  0x682e6ff3U,
  0x748f82eeU,
  0x78a5636fU,
  0x84c87814U,
  0x8cc70208U,
  0x90befffaU,
  0xa4506cebU, // 60
  0xbef9a3f7U,
  0xc67178f2U,
  0x98c7e2a2U,
  0xfc08884dU,
  0xcd2a11aeU,
  0x510e527fU,
  0x9b05688cU,
  0xC3910C8EU,
  0xfb6feee7U,
  0x2a01a605U, // 70
  0x0c2e12e0U,
  0x4498517BU,
  0x6a09e667U,
  0xa4ce148bU,
  0x95F61999U,
  0xc19bf174U,
  0xBB67AE85U,
  0x3C6EF372U,
  0xA54FF53AU,
  0x1F83D9ABU, // 80
  0x5BE0CD19U,
  0x5C5C5C5CU,
  0x36363636U,
  0x80000000U,
  0x000003FFU,
  0x00000280U,
  0x000004a0U,
  0x00000300U
};

#define rotl(x,y) rotate(x,y)
#define Ch(x,y,z) bitselect(z,y,x)
#define Maj(x,y,z) Ch((x^z),y,z)

#define EndianSwap(n) (rotl(n & ES[0], 24U)|rotl(n & ES[1], 8U))

#define Tr2(x)    (rotl(x, 30U) ^ rotl(x, 19U) ^ rotl(x, 10U))
#define Tr1(x)    (rotl(x, 26U) ^ rotl(x, 21U) ^ rotl(x, 7U))
#define Wr2(x)    (rotl(x, 25U) ^ rotl(x, 14U) ^ (x>>3U))
#define Wr1(x)    (rotl(x, 15U) ^ rotl(x, 13U) ^ (x>>10U))

#define RND(a, b, c, d, e, f, g, h, k)  \
h += Tr1(e);      \
h += Ch(e, f, g);     \
h += k;       \
d += h;       \
h += Tr2(a);      \
h += Maj(a, b, c);

void sha256(uint4*restrict state0,uint4*restrict state1, const uint4 block0, const uint4 block1, const uint4 block2, const uint4 block3)
{
  uint4 S0 = *state0;
  uint4 S1 = *state1;
  
  #define A S0.x
  #define B S0.y
  #define C S0.z
  #define D S0.w
  #define E S1.x
  #define F S1.y
  #define G S1.z
  #define H S1.w
  
  uint4 W[4];
  
  W[ 0].x = block0.x;
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[0]);
  W[ 0].y = block0.y;
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[1]);
  W[ 0].z = block0.z;
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[2]);
  W[ 0].w = block0.w;
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[3]);
  
  W[ 1].x = block1.x;
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[4]);
  W[ 1].y = block1.y;
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[5]);
  W[ 1].z = block1.z;
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[6]);
  W[ 1].w = block1.w;
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[7]);
  
  W[ 2].x = block2.x;
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[8]);
  W[ 2].y = block2.y;
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[9]);
  W[ 2].z = block2.z;
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[10]);
  W[ 2].w = block2.w;
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[11]);
  
  W[ 3].x = block3.x;
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[12]);
  W[ 3].y = block3.y;
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[13]);
  W[ 3].z = block3.z;
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[14]);
  W[ 3].w = block3.w;
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[76]);
  
  W[ 0].x += Wr1(W[ 3].z) + W[ 2].y + Wr2(W[ 0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[15]);
  
  W[ 0].y += Wr1(W[ 3].w) + W[ 2].z + Wr2(W[ 0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[16]);
  
  W[ 0].z += Wr1(W[ 0].x) + W[ 2].w + Wr2(W[ 0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[17]);
  
  W[ 0].w += Wr1(W[ 0].y) + W[ 3].x + Wr2(W[ 1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[18]);
  
  W[ 1].x += Wr1(W[ 0].z) + W[ 3].y + Wr2(W[ 1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[19]);
  
  W[ 1].y += Wr1(W[ 0].w) + W[ 3].z + Wr2(W[ 1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[20]);
  
  W[ 1].z += Wr1(W[ 1].x) + W[ 3].w + Wr2(W[ 1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[21]);
  
  W[ 1].w += Wr1(W[ 1].y) + W[ 0].x + Wr2(W[ 2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[22]);
  
  W[ 2].x += Wr1(W[ 1].z) + W[ 0].y + Wr2(W[ 2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[23]);
  
  W[ 2].y += Wr1(W[ 1].w) + W[ 0].z + Wr2(W[ 2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[24]);
  
  W[ 2].z += Wr1(W[ 2].x) + W[ 0].w + Wr2(W[ 2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[25]);
  
  W[ 2].w += Wr1(W[ 2].y) + W[ 1].x + Wr2(W[ 3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[26]);
  
  W[ 3].x += Wr1(W[ 2].z) + W[ 1].y + Wr2(W[ 3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[27]);
  
  W[ 3].y += Wr1(W[ 2].w) + W[ 1].z + Wr2(W[ 3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[28]);
  
  W[ 3].z += Wr1(W[ 3].x) + W[ 1].w + Wr2(W[ 3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[29]);
  
  W[ 3].w += Wr1(W[ 3].y) + W[ 2].x + Wr2(W[ 0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[30]);
  
  W[ 0].x += Wr1(W[ 3].z) + W[ 2].y + Wr2(W[ 0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[31]);
  
  W[ 0].y += Wr1(W[ 3].w) + W[ 2].z + Wr2(W[ 0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[32]);
  
  W[ 0].z += Wr1(W[ 0].x) + W[ 2].w + Wr2(W[ 0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[33]);
  
  W[ 0].w += Wr1(W[ 0].y) + W[ 3].x + Wr2(W[ 1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[34]);
  
  W[ 1].x += Wr1(W[ 0].z) + W[ 3].y + Wr2(W[ 1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[35]);
  
  W[ 1].y += Wr1(W[ 0].w) + W[ 3].z + Wr2(W[ 1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[36]);
  
  W[ 1].z += Wr1(W[ 1].x) + W[ 3].w + Wr2(W[ 1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[37]);
  
  W[ 1].w += Wr1(W[ 1].y) + W[ 0].x + Wr2(W[ 2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[38]);
  
  W[ 2].x += Wr1(W[ 1].z) + W[ 0].y + Wr2(W[ 2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[39]);
  
  W[ 2].y += Wr1(W[ 1].w) + W[ 0].z + Wr2(W[ 2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[40]);
  
  W[ 2].z += Wr1(W[ 2].x) + W[ 0].w + Wr2(W[ 2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[41]);
  
  W[ 2].w += Wr1(W[ 2].y) + W[ 1].x + Wr2(W[ 3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[42]);
  
  W[ 3].x += Wr1(W[ 2].z) + W[ 1].y + Wr2(W[ 3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[43]);
  
  W[ 3].y += Wr1(W[ 2].w) + W[ 1].z + Wr2(W[ 3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[44]);
  
  W[ 3].z += Wr1(W[ 3].x) + W[ 1].w + Wr2(W[ 3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[45]);
  
  W[ 3].w += Wr1(W[ 3].y) + W[ 2].x + Wr2(W[ 0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[46]);
  
  W[ 0].x += Wr1(W[ 3].z) + W[ 2].y + Wr2(W[ 0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[47]);
  
  W[ 0].y += Wr1(W[ 3].w) + W[ 2].z + Wr2(W[ 0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[48]);
  
  W[ 0].z += Wr1(W[ 0].x) + W[ 2].w + Wr2(W[ 0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[49]);
  
  W[ 0].w += Wr1(W[ 0].y) + W[ 3].x + Wr2(W[ 1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[50]);
  
  W[ 1].x += Wr1(W[ 0].z) + W[ 3].y + Wr2(W[ 1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[51]);
  
  W[ 1].y += Wr1(W[ 0].w) + W[ 3].z + Wr2(W[ 1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[52]);
  
  W[ 1].z += Wr1(W[ 1].x) + W[ 3].w + Wr2(W[ 1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[53]);
  
  W[ 1].w += Wr1(W[ 1].y) + W[ 0].x + Wr2(W[ 2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[54]);
  
  W[ 2].x += Wr1(W[ 1].z) + W[ 0].y + Wr2(W[ 2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[55]);
  
  W[ 2].y += Wr1(W[ 1].w) + W[ 0].z + Wr2(W[ 2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[56]);
  
  W[ 2].z += Wr1(W[ 2].x) + W[ 0].w + Wr2(W[ 2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[57]);
  
  W[ 2].w += Wr1(W[ 2].y) + W[ 1].x + Wr2(W[ 3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[58]);
  
  W[ 3].x += Wr1(W[ 2].z) + W[ 1].y + Wr2(W[ 3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[59]);
  
  W[ 3].y += Wr1(W[ 2].w) + W[ 1].z + Wr2(W[ 3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[60]);
  
  W[ 3].z += Wr1(W[ 3].x) + W[ 1].w + Wr2(W[ 3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[61]);
  
  W[ 3].w += Wr1(W[ 3].y) + W[ 2].x + Wr2(W[ 0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[62]);
  
  #undef A
  #undef B
  #undef C
  #undef D
  #undef E
  #undef F
  #undef G
  #undef H
  
  *state0 += S0;
  *state1 += S1;
}

void SHA256_fresh(uint4*restrict state0,uint4*restrict state1, const uint4 block0, const uint4 block1, const uint4 block2, const uint4 block3)
{
  #define A (*state0).x
  #define B (*state0).y
  #define C (*state0).z
  #define D (*state0).w
  #define E (*state1).x
  #define F (*state1).y
  #define G (*state1).z
  #define H (*state1).w
  
  uint4 W[4];
  
  W[0].x = block0.x;
  D= K[63] +W[0].x;
  H= K[64] +W[0].x;
  
  W[0].y = block0.y;
  C= K[65] +Tr1(D)+Ch(D, K[66], K[67])+W[0].y;
  G= K[68] +C+Tr2(H)+Ch(H, K[69] ,K[70]);
  
  W[0].z = block0.z;
  B= K[71] +Tr1(C)+Ch(C,D,K[66])+W[0].z;
  F= K[72] +B+Tr2(G)+Maj(G,H, K[73]);
  
  W[0].w = block0.w;
  A= K[74] +Tr1(B)+Ch(B,C,D)+W[0].w;
  E= K[75] +A+Tr2(F)+Maj(F,G,H);
  
  W[1].x = block1.x;
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[4]);
  W[1].y = block1.y;
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[5]);
  W[1].z = block1.z;
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[6]);
  W[1].w = block1.w;
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[7]);
  
  W[2].x = block2.x;
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[8]);
  W[2].y = block2.y;
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[9]);
  W[2].z = block2.z;
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[10]);
  W[2].w = block2.w;
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[11]);
  
  W[3].x = block3.x;
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[12]);
  W[3].y = block3.y;
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[13]);
  W[3].z = block3.z;
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[14]);
  W[3].w = block3.w;
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[76]);
  
  W[0].x += Wr1(W[3].z) + W[2].y + Wr2(W[0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[15]);
  
  W[0].y += Wr1(W[3].w) + W[2].z + Wr2(W[0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[16]);
  
  W[0].z += Wr1(W[0].x) + W[2].w + Wr2(W[0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[17]);
  
  W[0].w += Wr1(W[0].y) + W[3].x + Wr2(W[1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[18]);
  
  W[1].x += Wr1(W[0].z) + W[3].y + Wr2(W[1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[19]);
  
  W[1].y += Wr1(W[0].w) + W[3].z + Wr2(W[1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[20]);
  
  W[1].z += Wr1(W[1].x) + W[3].w + Wr2(W[1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[21]);
  
  W[1].w += Wr1(W[1].y) + W[0].x + Wr2(W[2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[22]);
  
  W[2].x += Wr1(W[1].z) + W[0].y + Wr2(W[2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[23]);
  
  W[2].y += Wr1(W[1].w) + W[0].z + Wr2(W[2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[24]);
  
  W[2].z += Wr1(W[2].x) + W[0].w + Wr2(W[2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[25]);
  
  W[2].w += Wr1(W[2].y) + W[1].x + Wr2(W[3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[26]);
  
  W[3].x += Wr1(W[2].z) + W[1].y + Wr2(W[3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[27]);
  
  W[3].y += Wr1(W[2].w) + W[1].z + Wr2(W[3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[28]);
  
  W[3].z += Wr1(W[3].x) + W[1].w + Wr2(W[3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[29]);
  
  W[3].w += Wr1(W[3].y) + W[2].x + Wr2(W[0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[30]);
  
  W[0].x += Wr1(W[3].z) + W[2].y + Wr2(W[0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[31]);
  
  W[0].y += Wr1(W[3].w) + W[2].z + Wr2(W[0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[32]);
  
  W[0].z += Wr1(W[0].x) + W[2].w + Wr2(W[0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[33]);
  
  W[0].w += Wr1(W[0].y) + W[3].x + Wr2(W[1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[34]);
  
  W[1].x += Wr1(W[0].z) + W[3].y + Wr2(W[1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[35]);
  
  W[1].y += Wr1(W[0].w) + W[3].z + Wr2(W[1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[36]);
  
  W[1].z += Wr1(W[1].x) + W[3].w + Wr2(W[1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[37]);
  
  W[1].w += Wr1(W[1].y) + W[0].x + Wr2(W[2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[38]);
  
  W[2].x += Wr1(W[1].z) + W[0].y + Wr2(W[2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[39]);
  
  W[2].y += Wr1(W[1].w) + W[0].z + Wr2(W[2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[40]);
  
  W[2].z += Wr1(W[2].x) + W[0].w + Wr2(W[2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[41]);
  
  W[2].w += Wr1(W[2].y) + W[1].x + Wr2(W[3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[42]);
  
  W[3].x += Wr1(W[2].z) + W[1].y + Wr2(W[3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[43]);
  
  W[3].y += Wr1(W[2].w) + W[1].z + Wr2(W[3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[44]);
  
  W[3].z += Wr1(W[3].x) + W[1].w + Wr2(W[3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[45]);
  
  W[3].w += Wr1(W[3].y) + W[2].x + Wr2(W[0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[46]);
  
  W[0].x += Wr1(W[3].z) + W[2].y + Wr2(W[0].y);
  RND(A,B,C,D,E,F,G,H, W[0].x+ K[47]);
  
  W[0].y += Wr1(W[3].w) + W[2].z + Wr2(W[0].z);
  RND(H,A,B,C,D,E,F,G, W[0].y+ K[48]);
  
  W[0].z += Wr1(W[0].x) + W[2].w + Wr2(W[0].w);
  RND(G,H,A,B,C,D,E,F, W[0].z+ K[49]);
  
  W[0].w += Wr1(W[0].y) + W[3].x + Wr2(W[1].x);
  RND(F,G,H,A,B,C,D,E, W[0].w+ K[50]);
  
  W[1].x += Wr1(W[0].z) + W[3].y + Wr2(W[1].y);
  RND(E,F,G,H,A,B,C,D, W[1].x+ K[51]);
  
  W[1].y += Wr1(W[0].w) + W[3].z + Wr2(W[1].z);
  RND(D,E,F,G,H,A,B,C, W[1].y+ K[52]);
  
  W[1].z += Wr1(W[1].x) + W[3].w + Wr2(W[1].w);
  RND(C,D,E,F,G,H,A,B, W[1].z+ K[53]);
  
  W[1].w += Wr1(W[1].y) + W[0].x + Wr2(W[2].x);
  RND(B,C,D,E,F,G,H,A, W[1].w+ K[54]);
  
  W[2].x += Wr1(W[1].z) + W[0].y + Wr2(W[2].y);
  RND(A,B,C,D,E,F,G,H, W[2].x+ K[55]);
  
  W[2].y += Wr1(W[1].w) + W[0].z + Wr2(W[2].z);
  RND(H,A,B,C,D,E,F,G, W[2].y+ K[56]);
  
  W[2].z += Wr1(W[2].x) + W[0].w + Wr2(W[2].w);
  RND(G,H,A,B,C,D,E,F, W[2].z+ K[57]);
  
  W[2].w += Wr1(W[2].y) + W[1].x + Wr2(W[3].x);
  RND(F,G,H,A,B,C,D,E, W[2].w+ K[58]);
  
  W[3].x += Wr1(W[2].z) + W[1].y + Wr2(W[3].y);
  RND(E,F,G,H,A,B,C,D, W[3].x+ K[59]);
  
  W[3].y += Wr1(W[2].w) + W[1].z + Wr2(W[3].z);
  RND(D,E,F,G,H,A,B,C, W[3].y+ K[60]);
  
  W[3].z += Wr1(W[3].x) + W[1].w + Wr2(W[3].w);
  RND(C,D,E,F,G,H,A,B, W[3].z+ K[61]);
  
  W[3].w += Wr1(W[3].y) + W[2].x + Wr2(W[0].x);
  RND(B,C,D,E,F,G,H,A, W[3].w+ K[62]);
  
  #undef A
  #undef B
  #undef C
  #undef D
  #undef E
  #undef F
  #undef G
  #undef H
  
  *state0 += (uint4)(K[73], K[77], K[78], K[79]);
  *state1 += (uint4)(K[66], K[67], K[80], K[81]);
}

void sha256SwapByteOrder(uint4 *data)
{
  *data = (uint4){
    EndianSwap((*data).x),
    EndianSwap((*data).y),
    EndianSwap((*data).z),
    EndianSwap((*data).w)
  };
}