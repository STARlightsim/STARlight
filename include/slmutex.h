#ifndef SLMUTEX_H
#define SLMUTEX_H

#include <pthread.h>

class MutexPosix
{
private:
 
  MutexPosix(const MutexPosix &);
  
  pthread_mutex_t _mutex;
  
public:
 
  MutexPosix() {}
  int lock() { return pthread_mutex_lock(&_mutex); }
  int unlock() { return pthread_mutex_unlock(&_mutex); }
  
};

template<typename M>
class Lockguard
{
private:
  
  M *_mutex;
  Lockguard(const Lockguard & guard); // Do not implement
 
public:
  
  Lockguard(M *mutex): _mutex(mutex) {_mutex->lock();}
  ~Lockguard() { _mutex->unlock(); }
  
};

#endif