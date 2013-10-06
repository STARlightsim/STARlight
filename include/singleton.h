#ifndef SINGLETON_H
#define SINGLETON_H

#ifdef CPP11
#include <atomic>
#include <mutex>


template<class T>
class Singleton {
public:
  static T& instance()
  {
    T * tmp = instance_.load(std::memory_order_consume);
    if (!tmp) {
      std::lock_guard<std::mutex> guard(instantiation_mutex);
      tmp = instance_.load(std::memory_order_consume);
      if (!tmp) {
        tmp = new T;
        instance_.store(tmp, std::memory_order_release);
      }
    }
    return *tmp;
  }
private:
  static std::atomic<T *> instance_;
  static std::mutex instantiation_mutex;
};

template<class T>
std::atomic<T *> Singleton<T>::instance_(0);

template<class T>
std::mutex Singleton<T>::instantiation_mutex;

#else

#include "slmutex.h"

template<class T>
class Singleton {
public:
  static T& instance()
  {
    T *tmp = _instance;
    if (!tmp) {
      Lockguard<MutexPosix> guard(&_instantiation_mutex);
      if (!tmp) {
        tmp = new T;
	_instance = tmp;
      }
    }
    return *tmp;
  }
private:
  
  static T * _instance;
  
  static MutexPosix _instantiation_mutex;
  
};

template<class T>
T *Singleton<T>::_instance(0);

template<class T>
MutexPosix Singleton<T>::_instantiation_mutex;

#endif

#endif

  