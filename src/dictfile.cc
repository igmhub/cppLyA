#include "dictfile.h"
#include "fileutils.h"
#include "fatalerror.h"
#include <stdio.h>
#include <fstream>

DictFileEntry::DictFileEntry(DictFile& F)
  : file(F)
{
  // maybe we should update the elements here.
}

DictFileEntry::DictFileEntry(const char *Line, DictFile &F)
  : file(F) 
{
  DecomposeString(elements, Line, " ");
}

DictFileEntry::Val DictFileEntry::Value(const string &Key, 
				       const bool DiesIfAbsent) const
{
  int offset = file.Dict().Locate(Key);
  if (offset >= 0 && offset < int(elements.size())) return Val(elements[offset]);
  if (DiesIfAbsent) FatalError("did not find Key "+Key+" in file "
			       +file.FileName(), true);
  cout << "return 'empty' string" << endl;
  return Val("empty");
}

DictFileEntry::Val DictFileEntry::Value(const unsigned Index, 
				       const bool DiesIfAbsent) const
{
  if (Index < elements.size()) return Val(elements[Index]);
  if (DiesIfAbsent) FatalError("did not find Key at required index in file "
			       +file.FileName());
  cout << "return 'empty' string" << endl;
  return Val("empty");
}



bool DictFileEntry::HasKey(const string &Key) const
{
  return file.HasKey(Key);
}

      
void DictFileEntry::AddKey(const string &Key, const string &Val)
{
  if (!file.HasKey(Key)) file.AddKey(Key);
  int offset = file.Dict().Locate(Key);
  if (offset == int(elements.size()))
    elements.push_back(Val);
  else
    cerr << " DicFileEntry::AddKey problem with Key " << Key 
	 << " Val " << Val << endl;
}

void DictFileEntry::AddKey(const string &Key, const double &Val)
{
  char s[80];
  sprintf(s,"%12.12E",Val);
  AddKey(Key,string(s));
}


void DictFileEntry::ModKey(const string &Key, const string &Val)
{
  if (!file.HasKey(Key)) {
    cerr << " DicFileEntry::ModKey no such key " << Key << " use AddKey" << endl;
    return;
  }
  int offset = file.Dict().Locate(Key);
  elements[offset]=Val;
}

void DictFileEntry::ModKey(const string &Key, const double &Val)
{
  char s[80];
  //sprintf(s,"%f",Val);
  sprintf(s,"%.10f",Val);
  ModKey(Key,string(s));
}



void DictFileEntry::writen(ofstream & pr) const
{
  for (int ii = 0 ; ii < file.Dict().size() ; ii++)
    pr << elements[ii] << " " ;

}



#include <cstring>

DictFile::DictFile()
  : fileName("")
{
  // no entries
}

void DictFile::Read(const string &FileName)
{
  fileName=  FileName ;
  ifstream f(FileName.c_str());
  if (!f) FatalError(" cannot open "+FileName);
 
  int Nmax = 10000 ;
  char * line = new char[Nmax];
  bool foundEnd = false;
  bool foundMat = false;
  while (f.getline(line,Nmax))
    {
      // skip leading spaces
      char *start = line;
      while (*start == ' ') start++;
      // cut at \n
      char *cr = rindex(start,'\n'); if (cr) *cr = '\0';

      if (*start == '\0') continue;
      if (*start == '#')
	{
	  start ++; // skip '#'
	  while (*start == ' ') start++;
	  if (strstr(start,"end") == start) {foundEnd = true; continue;}
	  if (strstr(start,"MATRIX") == start) {foundMat = true; break;}
	  char *column = strchr(start,':');
	  if (column)
	    {
	      if (foundEnd)
		{
		  std::cerr << " bizarre file structure : found tags after end"
			    << std::endl;
		}
	      column = start + strcspn(start," \t:");
	      *column = '\0';
	      int presentSize = dict.size();
	      string tag(start);
	      dict[tag] = presentSize;
	      // DEBUG
	      //	      cout << "size " << presentSize << ' ' << tag << endl;
	    }
	}
      else if (*start == '@') 
	{
	  GlobalKeys.ProcessLine(line);
	} 
      else 
	{ 
	  push_back(DictFileEntry(start,*this));
	}
    }
  if(foundMat) mat.readASCII(f);
  f.close();
  delete [] line ;
}

DictFile::DictFile(const string &FileName)
{
  Read(FileName);
}
void DictFile::Append(const string &FileName) 
{
  FILE *f = fopen(FileName.c_str(),"r");
  if (!f) FatalError(" cannot open "+FileName);
 
  int Nmax = 10000 ;
  char * line = new char[Nmax];
  bool foundEnd = false;
  while (fgets(line,Nmax,f))
    {
      // skip leading spaces
      char *start = line;
      while (*start == ' ') start++;
      // cut at \n
      char *cr = rindex(start,'\n'); if (cr) *cr = '\0';

      if (*start == '\0') continue;
      if (*start == '#')continue;
      if (*start == '@')continue;
      push_back(DictFileEntry(start,*this));
    }
  fclose(f);
  delete [] line ;
}

void DictFile::DumpKeys() const 
{
  for(Dictionnary::const_iterator it = dict.begin(); it!=dict.end(); ++it) {
    cout << it->second << " " << it->first << endl;
  }
}


bool DictFile::WriteHeader_(ofstream & pr, string suffixe) const {

  unsigned presentSize = dict.size();
 // invert the dictionnary:
  map<int,string> tags;
  for (Dictionnary::const_iterator it = dict.begin(); it != dict.end(); ++it)
    tags[it->second] = it->first;

  //write the header
  for (unsigned i = 0; i < presentSize; ++i)
    pr << "# " <<  tags[i] << suffixe << " : " << endl ;
  return(true);
}

bool DictFile::Write(const std::string &FileName, bool name_at_the_end) const

{
  size_t presentSize = dict.size();
  FILE * f = fopen(FileName.c_str(), "w");
  if (!f)
    {
      cerr << " DistFile::Write() : could not open " << FileName << endl;
      return false;
    }

  // write global keys and values
  vector<string> v(GlobalKeys.OutputLines());
  for (size_t k=0; k<v.size(); ++k) 
    fprintf(f,"@%s\n",v[k].c_str());
  
  // invert the dictionnary:
  map<int,string> tags;
  for (Dictionnary::const_iterator it = dict.begin(); it != dict.end(); ++it)
    tags[it->second] = it->first;

  //write the header
  bool name_seen = false ;
  for (unsigned i = 0; i < presentSize; ++i)
    {
      if (! name_at_the_end)
	fprintf(f,"# %s :\n",tags[i].c_str());
      else
	{
	  if(tags[i] != "name") 
	    fprintf(f,"# %s :\n",tags[i].c_str());
	  else
	    name_seen = true ;
	}
    }
  if (name_at_the_end && name_seen )
    fprintf(f,"# name :\n");
  fprintf(f,"# end\n");
  // write data
  for (const_iterator it = begin(); it != end(); ++it)
    {
      const DictFileEntry &entry = *it;
      for (size_t i =0; i < presentSize; ++i)
	{
	  const string& val(entry.Value(tags[i]));
	  if ((! name_at_the_end) || (tags[i] != "name")) 
	    fprintf(f,"%s ", val.c_str());
	}
      if (name_at_the_end && name_seen )
	{
	  const string& val(entry.Value("name"));
	  fprintf(f,"%s ", val.c_str());
	}
      fprintf(f,"\n");
    }
  fclose(f);
  return true;
}

void DictFile::AddKey(const string &Key)
{
  if (!HasKey(Key))
    {
      int presentSize = dict.size();
      dict[Key] = presentSize;
    }
}




string DictFile::GlobalValue(const string &Key, const bool ExitIfAbsent) const
{			  
  if(! GlobalKeys.HasKey(Key)) {
    if (ExitIfAbsent) FatalError(" giving up : we miss key "+Key+" in "+fileName);
    return "";
  }
  return GlobalKeys.getStringValue(Key);
}


void DictFile::RemoveGlobalKey( const string &Key) {
  if( HasGlobalKey(Key) ) {
    GlobalKeys.Erase(Key);
  }else{
    cerr << "warning no such key " << Key << " in global values of " << fileName << endl;
  }
}

void DictFile::AddGlobalKey( const string &Key, const string &Val) {
  if( HasGlobalKey(Key) ) {
    cerr << "warning overwriting " << Key << "=" << GlobalValue(Key) << " in " << fileName << endl;
  }
  GlobalKeys.AddKey(Key,Val);
}


std::vector<std::string> DictFile::Keys() const {
  std::vector<std::string> keys;
  for(std::map<std::string,int>::const_iterator it=dict.begin();it!=dict.end();++it)
    keys.push_back(it->first);
  return keys;
}
