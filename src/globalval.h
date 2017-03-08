// -*- C++ -*-
// 
// \file globalval.h
// 
// Last modified: $Date: 2010-10-07 15:35:16 $
// By:            $Author: guy $
// 

#ifndef GLOBALVAL__H
#define GLOBALVAL__H


#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

//! to store in files things like "Key value(s)" things.
class GlobalVal : public  map<string, vector<string> > {

private :
  
  //! default=false , true when a key @ALLOW_MULTIPLE_TAGS is read
  bool allow_multiple_tags;

public :

  GlobalVal() : allow_multiple_tags(false) {};

  bool AddKey(const string &Key, const vector<string> &Values);
  
  bool AddKey(const string &Key, const string &Value);
  
  bool AddKey(const string& Key, const list<string>& Values);
  
  bool AddKey(const string &Key, const vector<double> &Values);
  
  bool AddKey(const string &Key, const double &Value);

  bool AddKey(const string& Key, const list<double>& Values);
  
  size_t NKey() const;

  bool HasKey(const string &Key) const;
  
  string         getStringValue(const string& Key, bool fatal_if_failed = true) const;
  
  vector<string> getStringValues(const string& Key, bool fatal_if_failed = true) const;
  
  double         getDoubleValue(const string& Key, bool fatal_if_failed = true) const;
  
  vector<double> getDoubleValues(const string& Key, bool fatal_if_failed = true) const;

  vector<string> OutputLines() const;

  vector<string> Keys() const;

  void Erase(const string& Key);
  
  bool ProcessLine(char *Line);
  
  void Read(const std::string &FileName, bool break_when_digit=true);
  // this constructor is to be used when you only need to read '@' lines in a file
  GlobalVal(const std::string &FileName, bool break_when_digit=true);

  //  GlobalVal(const GlobalVal &O) { cout << " on copie un globalval " << endl;}


private:
  template<typename T>
  bool GenericAddKey(const string& Key, const T& Values);
};

std::ostream& operator <<(ostream &S, const GlobalVal &G);

#endif /* GLOBALVAL__H */
