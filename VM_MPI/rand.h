#ifndef RAND_H
#define RAND_H

// Generate a random number
struct Ran
{
  unsigned long long u, v, w;
  Ran(unsigned long long j) :v(4101842887655102017LL), w(1) {
    u = j^v; int64();
    v = u; int64();
    w = v; int64();
  }
  inline unsigned long long int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }
  inline double doub() {
    return 5.42101086242752217E-20 * int64();
  }
  inline unsigned int int32() {
    return (unsigned int)int64();
  }
};


template<class T>
void shuffle(T *a, int n, Ran *myran) {
  for (int i = n - 1; i >= 0; i--) {
    // generate a random int j that 0 <= j <= i  
    int j = int(myran->doub() * (i + 1));
    if (j > i)
      j = i;
    else if (j < 0)
      j = 0;
    T tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }
}
#endif

