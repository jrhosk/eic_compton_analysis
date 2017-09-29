#include <iostream>
#include "MsgStream.hh"

namespace Sys {
  extern MsgStream::SystemMessage SysMsg;
  extern MsgStream::SystemError SysError;

  static std::ostream& endl(std::ostream& stream)
  {
    stream << kWhiteFont << std::endl;
    
    return stream;
  }
}
