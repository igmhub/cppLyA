#include <iostream>
#include "pathutils.h"
#include "fileutils.h"
#include "datacards.h"

using namespace std;

string AssembleFileName(const string &DirName, const string &FileName,
			const bool Warn)
{
  ///// if (FileExists(FileName)) return FileName;
  string tmp = AddSlash(DirName);
  string s1 = tmp+FileName;
  if (FileExists(s1)) { 
    return s1;
  }
  string s2 = tmp+BaseName(FileName);
  if (Warn && !FileExists(s2))
    {
      cerr << " cannot find '" << FileName << "' nor '" << s1 << "' nor '" << s2 << "'" << endl;
      cerr <<" expect troubles" << endl;
    }
  return s2;
}

string GetDataPath() {
  static char* salt_path = NULL;
  if(salt_path) return salt_path;
  
  
  salt_path = getenv("SALTPATH");
  
  // get PATHMODEL anyway to write a warning message if exists and give different path 
  char* path_model = getenv("PATHMODEL"); 
  
  if(path_model) {
    
    if(salt_path) {
      if(AddSlash(path_model) != AddSlash(salt_path)) {
	cerr << "ERROR conflicting (new) SALTPATH='" << salt_path << "' and (old) PATHMODEL='" << path_model << "'" << endl;
	cerr << "Please unset one of the two (SALTPATH is more recent environment variable)" << endl;
	cerr << "Exiting now to avoid further trouble" << endl;
	exit(12);
      }
    } else { /* not salt_path */
      salt_path = path_model;
      cerr << "WARNING  please use (new) environment variable SALTPATH instead of PATHMODEL='" << path_model << "' in the future." << endl;
    }
  }
  
  if(salt_path) return salt_path;
  
  cerr << " you should setenv SALTPATH (new name for environment PATHMODEL (which is not set either)) to the location of calib files" << endl;
  
  exit(12);
  

}


string LocateFile(const string &FileName)
{
  return AssembleFileName(GetDataPath(),FileName);
}

string LocateFileFromCard(const string &Key)
{
  string dirname = GetDataPath();
  string mainfile = dirname+"/fitmodel.card";
  DataCards cards(mainfile);
  if(! cards.HasKey(Key))
    std::cerr << " cannot find key " << Key << " in main file " << mainfile << endl;
  else
    {
      string FileName = cards.SParam(Key);
#ifdef DEBUG
      cout << "LocateFileFromCard: " << Key << " -> " << FileName << endl;
#endif
      return AssembleFileName(dirname,FileName);
    }
  return AssembleFileName(dirname,Key);
}



