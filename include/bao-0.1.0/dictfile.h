// -*- C++ -*-
// 
#ifndef DICTFILE__H
#define DICTFILE__H


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <globalval.h>

#include <cstdlib> // for atof

#include "matvect.h"

//using namespace std;

class DictFile;

class Dictionnary : public std::map<std::string,int> 
{ // just add a const [] to map
  public :
  int Locate(const std::string &Key) const 
    {
      const_iterator i = find(Key);
      return (i == end()) ? -1 : i->second;
    } 
  std::map<int,std::string> Key_List()const {
        std::map<int,std::string> tags;
        for (Dictionnary::const_iterator it = begin(); it != end(); ++it)
            tags[it->second] = it->first;
        return tags;
    }

};

class DictFile;

//! only its Value(string Key) routine is useful
class DictFileEntry
{

private :
  std::vector<std::string> elements;
  DictFile &file;

public :
  DictFileEntry(DictFile& F);
  DictFileEntry(const char *Line, DictFile &F);


  struct Val
    {
      const std::string &val;
      Val( const std::string &Sval) : val(Sval) {};
      operator double() const { return atof(val.c_str());};
      operator std::string() const { return val;};
  };

  //! this routine enable double toto = line.Value("STUFF");
  Val Value(const std::string &Key, const bool DiesIfAbsent = true) const;
  Val Value(const unsigned Index, const bool DiesIfAbsent = true) const;
  bool HasKey(const std::string &Key) const;



  void AddKey(const std::string &Key, const std::string &Val);
  void AddKey(const std::string &Key, const double &Val);
  void ModKey(const std::string &Key, const std::string &Val);
  void ModKey(const std::string &Key, const double &Val);
  void writen(std::ofstream & pr) const;
  size_t size() const { return elements.size();}
  
};

  



/*! a class to access files like
# Flux : B Flux
# Fluxerr : error on B Flux
# FluxPsf : B FluxPsf
# FluxPsferr : error on B FluxPsf
# Day : Day
# Dayerr : error on Day
# AirMass : Airmass when relevant
# Absorption : Absorption when relevant
# Band : Band of obs
# Instrument : Instrument of abs
# format LightCurvePoint 1
# end
@BAND B
@INSTRUMENT toto
-13.5392 0.0509903 0 1e+30 -3.60747 0.01 0 1 B STANDARD
-13.5291 0.0412349 0 1e+30 -3.47747 0.01 0 1 B STANDARD
-13.4792 0.0412312 0 1e+30 -2.55747 0.01 0 1 B STANDARD
-13.4792 0.0509836 0 1e+30 -1.67747 0.01 0 1 B STANDARD
-13.4292 0.0510007 0 1e+30 0.442531 0.01 0 1 B STANDARD
*/
class DictFile : public std::list <DictFileEntry>
{
  Dictionnary dict;
  std::string fileName;
  Mat mat;

public :
  DictFile();
  DictFile(const std::string &FileName);
  void Read(const string &FileName);
  void Append(const std::string &FileName);
  const Dictionnary &Dict() const { return dict;}
  const std::string &FileName() const {return fileName;}
  bool HasKey(const std::string &Key) const { return dict.Locate(Key) != -1;}
  void AddKey(const std::string &Key);
  void RmKey(const std::string &Key);
  std::vector<std::string> Keys() const;
  void DumpKeys() const;
  
  GlobalVal GlobalKeys;

  bool HasGlobalKey( const std::string &Key) const 
  { return  GlobalKeys.HasKey(Key);};
  
  bool HasMat() const
  { return mat.SizeX() != 0;};
  
  Mat GetMat() const
  { return mat;};
  
  std::string GlobalValue(const std::string &Key, const bool FatalIfAbsent = true) const;
  void RemoveGlobalKey( const std::string &Key);
  void AddGlobalKey( const std::string &Key, const std::string &Val);
  
  bool WriteHeader_(ofstream & pr, string suffixe) const ;
  bool Write(const std::string &FileName, bool name_at_the_end = false) const;

};


typedef DictFile::const_iterator DictFileCIterator;
typedef DictFile::iterator DictFileIterator;


#ifdef EXAMPLE
example of usage of classe in this file :

DictList file(FileName);
double redshift = file.FlobalValue("REDSHIFT"); // dies if absent.
for (DictFileCIterator line = file.begin(); line != file.end(); 
  ++line)
{
  MyStruct a;
  if (line->HasKey("Stuff")) a->stuff = line->Value("Stuff");
}

#endif
    




#endif /* DICTFILE__H */
