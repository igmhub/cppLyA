// 
// \file globalval.cc
// 
// Last modified: $Date: 2010-10-07 15:35:16 $
// By:            $Author: guy $
// 
#include <iostream>
#include <sstream>
#include "globalval.h"
#include "fileutils.h" // for DecomposeString 

#include <cstdlib> // for atof
#include <cstring> // for strlen

bool GlobalVal::HasKey(const string &Key) const
{
  return (find(Key) != end());
}


size_t GlobalVal::NKey() const
{
  return this->size();
}


template<typename T>
bool GlobalVal::GenericAddKey(const string& Key, const T& Values) {
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
    

  typename T::const_iterator I;
  for(I=Values.begin();I!=Values.end();I++) {
    stringstream sstrm;    
    sstrm << *I;
    (*this)[Key].push_back(sstrm.str());
  }
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<string>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const vector<string> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
  
  (*this)[Key] = Values;
  return true;
}

bool GlobalVal::AddKey(const string &Key, const string &Value)
{
  vector<string> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}

void GlobalVal::Erase(const string &Key)
{
  map<string, vector<string> >::iterator it = find(Key);
  if(it!= end())
    erase(it);
}


bool GlobalVal::AddKey(const string &Key, const vector<double> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
  
  size_t i;
  for(i=0;i<Values.size();i++) {
    stringstream sstrm;
    sstrm << Values[i];
    (*this)[Key].push_back(sstrm.str());
  }
  //  (*this)[Key] = Values;
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<double>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const double &Value)
{
  vector<double> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}

static void warning_or_die(const string& message, bool fatal) {
  if(fatal)
    cerr << "FATAL ERROR";
  else
    cerr << "WARNING";
  cerr << message << endl;
  if(fatal)
    exit(12);
}

string GlobalVal::getStringValue(const string& Key, bool fatal_if_failed) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
       warning_or_die(" in GlobalVal::getStringValue, no key '"+Key+"'",fatal_if_failed);
      return "";
    }
  else return i->second[0];
}


vector<string> GlobalVal::getStringValues(const string& Key, bool fatal_if_failed) const
{
  const_iterator i = find(Key);
  if(i == end())
    {
      warning_or_die(" in GlobalVal::getStringValues, no key '"+Key+"'",fatal_if_failed);
      vector<string> ret;
      return ret;
    }
  return i->second;
}



double GlobalVal::getDoubleValue(const string &Key, bool fatal_if_failed) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      warning_or_die(" in GlobalVal::getDoubleValue, no key '"+Key+"'",fatal_if_failed);
      return 99;
    }
  else {
    char *endptr;
    double res = strtod(i->second[0].c_str(),&endptr);
    if(i->second[0].c_str() == endptr) {
      warning_or_die(" in GlobalVal::getDoubleValue, could not convert value '"
		     +i->second[0]+"' into double (for key '"+Key+"')",fatal_if_failed);
      return 99;
    }
   
    return res;

  }
}



vector<double> GlobalVal::getDoubleValues(const string &Key, bool fatal_if_failed) const
{
  vector<double> ret;
  const_iterator i = find(Key);
  if (i == end())
    {
      warning_or_die(" in GlobalVal::getDoubleValue, no key '"+Key+"'",fatal_if_failed);
      return ret;
    }
  size_t k;
  char *endptr;
  for(k=0;k<i->second.size();k++) {
    double res = strtod(i->second[k].c_str(),&endptr);
    if(i->second[k].c_str() == endptr) {
      warning_or_die(" in GlobalVal::getDoubleValue, could not convert value '"
		     +i->second[k]+"' into double (for key '"+Key+"')",fatal_if_failed);
      res = 99;
    } 
    ret.push_back(res);
  }
  return ret;
}


#include <sstream>

vector<string> GlobalVal::OutputLines() const
{
  vector<string> out;
  for (const_iterator i = begin(); i != end(); ++i)
    {
      //      const vector<double> &values = i->second;
      const vector<string> &values = i->second;
      if (values.size() == 0) continue;
      ostringstream s;
      s << i->first;
      for (size_t k=0; k < values.size(); ++k) s << ' ' << values[k];
      out.push_back(s.str());
    }
  return out;
}

vector<string> GlobalVal::Keys() const
{
  vector<string> out;
  for (const_iterator i = begin(); i != end(); ++i)
    out.push_back(i->first);
  return out;
}

std::ostream& operator <<(ostream &S, const GlobalVal &G)
{
  vector<string> v(G.OutputLines());
  for (size_t k=0; k<v.size(); ++k) S << "@" << v[k] << endl;
  // S << " glob size " << v.size() << endl;
  return S;
}
		   


#include <fstream>

void GlobalVal::Read(const std::string &FileName, bool break_when_digit) 
{
  
  FILE *f = fopen(FileName.c_str(),"r");
  if (!f) 
    {
      cerr << " GlobalVal : could not open " << FileName << endl;
      return;
    }
  char line[4096];
  while (fgets(line,4096,f)!=0)
    {
      if (line[0] == '@') 
	{
	  ProcessLine(line);
	  continue;
	}
      if (break_when_digit && isdigit(line[0])) break;
    }
  fclose(f);
}

GlobalVal::GlobalVal(const std::string &FileName, bool break_when_digit) : 
  allow_multiple_tags(false)
{
  Read(FileName,break_when_digit);
}


#include <stdlib.h>


bool GlobalVal::ProcessLine(char *Line)
{
  // use standard C stuff because it enables error checking
  size_t l = strlen(Line);
  if(l>0 && Line[l-1]=='\n') Line[l-1]='\0'; // remove '\n' at the end

  const char *s = Line;
  if (*s == '@') s++;

  char buf[256];
  if (sscanf(s,"%s",buf)!= 1)
    {
      cerr << " could no read a key in " << Line << endl;
      return false;
    }
  s += strlen(buf);
  string key(buf);

  if(key == "ALLOW_MULTIPLE_TAGS") {
    allow_multiple_tags = true;
    return true;
  }
  if (HasKey(key) && (!allow_multiple_tags))
    {
      cerr <<" already have a key labelled " << key << endl;
      return false;
    }

  vector<string> &values = (*this)[key];
  string strtmp = s;
  DecomposeString(values, strtmp, " ");
  return true;
}

