#ifndef message_hh
#define message_hh

#include <iostream>
#include <string>

#include "FontColor.hh"

class MsgStream : public std::ostream {

public:
  MsgStream();
  ~MsgStream();
  
  friend std::ostream& operator<<(std::ostream& os, MsgStream& t);
  
  struct SystemMessage
  {

    template <typename T>  SystemMessage& operator<<(const T& x)
    {
      std::cout << kBlueFont  << x;
      
      return *this;
    }
  };
  struct SystemError
  {

    template <typename T>  SystemError& operator<<(const T& x)
    {
      std::cerr << kRedFont  << x;
      
      return *this;
    }
  };
};  
#endif



