#include <stdint.h>
#include <stdio.h>

static const uint16_t rsqrt_table[32] = {
    65536, 46341, 32768, 23170, 16384,  /* 2^0 to 2^4 */
    11585,  8192,  5793,  4096,  2896,  /* 2^5 to 2^9 */
     2048,  1448,  1024,   724,   512,  /* 2^10 to 2^14 */
      362,   256,   181,   128,    90,  /* 2^15 to 2^19 */
       64,    45,    32,    23,    16,  /* 2^20 to 2^24 */
       11,     8,     6,     4,     3,  /* 2^25 to 2^29 */
        2,     1                         /* 2^30, 2^31 */
};

static int clz(uint32_t x){
  if(!x) return 32;
  int n = 0;
  if(!(x & 0xFFFF0000)) {n += 16; x <<= 16;}
  if(!(x & 0xFF000000)) {n += 8; x <<= 8;}
  if(!(x & 0xF0000000)) {n += 4; x <<= 4;}
  if(!(x & 0xC0000000)) {n += 2; x <<= 2;}
  if(!(x & 0x80000000)) {n += 1;}
  return n;
}

static uint64_t mul32(uint32_t a, uint32_t b){
  uint64_t r = 0;
  for(int i=0; i < 32; i++){
    if(b & (1U << i))
      r+=(uint64_t)a << i;
  }
  return r;
}


uint32_t fast_rsqrt(uint32_t x){
  if(x==0) return 0xFFFFFFFF;
  if(x==1) return 65536;
  
  int exp = 31-clz(x);
  
  uint32_t y = rsqrt_table[exp];
  
  if(x > (1u << exp)){
    uint32_t y_next = (exp < 31) ? rsqrt_table[exp + 1] : 0;
    uint32_t delta = y - y_next;
    uint32_t frac = (uint32_t)((((uint64_t)x - (1UL << exp))<< 16)>> exp);
    y-= (uint32_t)((delta * frac) >> 16);
  }
  
  for(int iter = 0; iter <2; iter++){
    uint32_t y2 = (uint32_t)mul32(y,y);
    uint32_t xy2 = (uint32_t)(mul32(x, y2) >> 16);
    y = (uint32_t)(mul32(y, (3u << 16) - xy2) >> 17);
  }
  
  return y;
}

uint32_t fast_distance_3d(int32_t dx, int32_t dy, int32_t dz)
{
    uint64_t dist_sq = (uint64_t)dx * dx +
                       (uint64_t)dy * dy +
                       (uint64_t)dz * dz;

    if (dist_sq > 0xFFFFFFFF)
        dist_sq >>= 16;

    uint32_t inv_dist = fast_rsqrt((uint32_t)dist_sq);

    /* sqrt(x) = x / sqrt(x) = x * (1/sqrt(x)) */
    uint64_t dist = ((uint64_t)dist_sq * inv_dist) >> 16;

    return (uint32_t)dist;
}

int main() {
    uint32_t x =UINT32_MAX;
    printf("Test1: %u\n", mul32(3,5));
    //printf("Test2: %u\n", fast_distance_3d(100, 200, 300));
    //printf("Test3: %u\n", fast_distance_3d(10000, 20000, 10000));
    return 0;
}