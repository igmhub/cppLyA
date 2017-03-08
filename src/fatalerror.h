#ifndef FATALERROR__H
#define FATALERROR__H

#include <string>

//! prints the message then exit() or abort() depending on Abort.
void FatalError(const std::string &Message, const bool Abort=false);


#endif
