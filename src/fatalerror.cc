#include "fatalerror.h"

#include <iostream> // for cerr
#include <stdlib.h> // for exit

void FatalError(const std::string &Message, const bool Abort)
{
  std::cerr << Message << std::endl << " we stop here ... " << std::endl;
  if (Abort) abort();
  exit(EXIT_FAILURE);
}
