#ifndef SNFITEXCEPTIONS_H
#define SNFITEXCEPTIONS_H

using namespace std;
#include <string>
#include <sstream>
#include <iostream>
#include <exception>
/* Exception base class

*/

/* Snfit exceptions rules (proposal?):

1) all exception classes should derive from SnfitException.
2) the message that this class handles should be designed as a "one liner"
   to be grep'ed in log files.
3) if you need to provide more information than a single line, just print 
   it out on std::cerr before throwing.
4) there is no built-in mechanism that ensures that the message is printed
   --> all mains should catch (SnfitException's) and PrintMessage(std::cerr)
   in order to make sure that the message is printed and the 
   failure status forwarded to the calling scripts.
5) if you catch, either you print or you (re) throw.


*/
 




//namespace toads {
#ifndef BuildExcMsg
#define BuildExcMsg(msg) Snfit::buildMessage(msg, __FILE__, __LINE__)
#endif


class SnfitException {
public:
  //! constructor
  SnfitException(string const& msg) : msg_(msg) {} 
  
  //! destructor -- does nothing. this class may (should) be derived.
  virtual ~SnfitException() {}
  
  //! append stuff to the initial message
  void   append(string const& msg) { msg_+=msg; }
  
  //! return the message
  //  string message() const { return msg_; }

  //! return the message
  const string& message() const { return msg_; }

  //! return the message, standard way
  virtual const char* what() const throw()
  {
    return msg_.c_str();
  }

  //! print message (here to force the "tag" string) for "grep" in logs
  void PrintMessage(ostream &s) const
  { s << "SNFIT_EXCEPTION: " << msg_ << endl; }

  
  //! uniform message format
  static string buildMessage(string msg, char const* file, size_t line) {
    stringstream ret;
    ret << "[" << file << "]" << "{" << line << "} " << msg;
    return ret.str();
  }

private:
  string msg_;
};



#endif


