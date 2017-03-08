#include <cstdio>
#include <cstdlib>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "fileutils.h"

#include <string>
#include <cctype>
#include <iostream>
#include <cstring> // for strcmp, strncmp ,strstr, ...


void DirName(const char *FullFileName, char *DirName)  /* DirName is assumed to point to an area at least as long as FullFileName */
{
  char *end;
  strcpy(DirName,FullFileName);
  end = DirName+strlen(DirName) - 1;
  if (*end == '/') *end = '\0'; /* ignore a trailing / as the dirname command */
  end = strrchr(DirName,'/');
  if (end) *end = '\0';
  else strcpy(DirName,".");
}

string DirName(const string &FullFileName)
{
  char toto[256];
  DirName(FullFileName.c_str(), toto);
  return string(toto);
}

void BaseName(const char *FullFileName, char *Name)
{
  char a_path[256];
  const char *begin;
  strcpy(a_path,FullFileName);
  if (a_path[strlen(a_path)-1] == '/') a_path[strlen(a_path)-1] = '\0';
  begin = strrchr(a_path,'/');
  if (!begin) strcpy(Name, a_path);
  else strcpy(Name,begin + 1);
}

string BaseName(const string &Name)
{
  char toto[512];
  BaseName(Name.c_str() ,toto);
  return string(toto);
}


string CutExtension(const string &Name)
{
  char result[512];
  strcpy(result,Name.c_str());
  char *dot = strrchr(result, '.');
  if (dot) *dot = '\0';
  return string(result);
}


int MKDir(const char *path, const bool Warn)
{
  if (FileExists(path)) return 1;
  char command[512];
  sprintf(command,"mkdir -p %s",path);
  int status  = system(command);
  // cout <<" debug :status for " << command << ' ' << status << endl;
  if (status && Warn)
    {
      cerr << " MKDir : could not create " << path << endl;
    }
  return (!status);
}
int cp(const char *path1, const char *path2, const bool Warn)
{
  char command[512];
  sprintf(command,"cp %s %s",path1,path2);
  int status  = system(command);
  // cout <<" debug :status for " << command << ' ' << status << endl;
  if (status && Warn)
    {
      cerr << " could not do " << command << endl;
    }
  return (!status);
}


string FullFileName(const string &FileName)
{
  if (strlen(FileName.c_str()) && FileName[0] == '/') return FileName;
  //char cwd[512];
  //getcwd(cwd,512); /* very nice : may not be equal to $cwd from shell (if there are links involved) */
  return  AddSlash(getenv("PWD")) + FileName;
}

static string RemoveDoubleSlash(const string &FileName)
{
  std::string result = FileName;
  // string::find does not tell clearly when the pattern is not found
  while (strstr(result.c_str(),"//"))
    {
      int pos = result.find("//");
      result.erase(pos,1);
    }
  return result;
}


string RelativeLink(const char *RealFile, const char *LinkLocation)
//  caveats : does not check if there are symbolic links in the given pathes.. 
//            misses physically identical pathes if they do not have the same name
//  this kind of utility routine should be stolen in shell or unix source code.
{
  string realFile = FullFileName(RealFile);
  string link     = FullFileName(LinkLocation);
  //
  realFile = RemoveDoubleSlash(realFile);
  link     = RemoveDoubleSlash(link);

 
  /* find the common part */
  int bound = min(realFile.length(), link.length());
  int i;
  for (i = 0; i< bound; ++i) 
    if (link[i] != realFile[i]) break;
  while (link[i] != '/') i--;
  i++;
  if (i == 1) /* they only share the leading '/' */ return realFile;
  /* assemble the result */  
  string result;
  for (const char* p = link.c_str() + i; (p=strchr(p,'/')) ; )
    { result += "../"; ++p; }
  result = result + (realFile.c_str()+i);
  return result;
}

int MakeRelativeLink(const std::string &RealFile, 
		     const std::string &LinkLocation)
{
  string link_value = RelativeLink(StandardPath(RealFile).c_str(), 
				   StandardPath(LinkLocation).c_str());
  string command = " ln -fs " + link_value + " " + string(LinkLocation);
  if (system(command.c_str()) == 0) return 1;
  cerr << " could not issue " << command << endl;
  return 0;
}

string AddSlash(const string &path)
{
  if (path[path.length()-1] == '/') return path ; else return path+"/";
}


int FileExists(const char *FileName)
{
  return (access(FileName,F_OK) == 0);
}

int FileExists(const string &FileName)
{
  return (access(FileName.c_str(),F_OK) == 0);
}

int FileIsWritable(const char *FileName)
{
  struct stat a_stat;
  if (stat(FileName, &a_stat) == 0) /* success ! */
    {
      return (a_stat.st_mode & (S_IWUSR|S_IWGRP|S_IWOTH));
    }
  return 0;
}

/* borrowed in ximtool/fitsio.c */
int  IsFits (const string &Name)
{
  FILE *fp;
  int value = 0;
  char keyw[8], val;

  if ((fp = fopen (Name.c_str(), "r"))) {
    fscanf (fp, "%6s = %c", keyw, &val);
    if (strcmp ("SIMPLE", keyw) == 0 && val == 'T')
      value = 1;
    fclose (fp);
  }
  return value;
}



#include <sys/types.h>
#include <dirent.h>



/* open a directory and returns the file names it contains */
int DirectoryContents(const string& DirName, StringList &FileNames)
{
  DIR *this_dir = opendir(DirName.c_str());
  FileNames.clear();
  if (!this_dir) return 0;
  struct dirent* entry;
  while ((entry=readdir(this_dir)))
    {
      if (!entry->d_name) continue;
      if (strcmp(entry->d_name,".") == 0 || strcmp(entry->d_name,"..") == 0) continue;
      FileNames.push_back(entry->d_name);
    }
  FileNames.sort();
  return 1;
}


#ifdef I_DO_NOT_KNOW_HOW_TO_DO
#include <sys/stat.h>
#include <unistd.h>

int IsLink(const string &FileName)
{
  struct stat stat_struct;
  /* lstat return 0 on success */
  return (lstat(FileName.c_str(), &stat_struct) ==0  && 
	  S_ISLNK(stat_struct.st_mode));

}

#endif /* I_DO_NOT_KNOW_HOW_TO_DO */


int IsDirectory(const string &FileName)
{
  struct stat stat_struct;
  /* stat return 0 on success */
  if (stat(FileName.c_str(), &stat_struct) != 0) return 0;
  return S_ISDIR(stat_struct.st_mode);
}

//! return true for filename ending by ".Z" or ".gz"
bool IsCompressed(const string &FileName)
{
  const char*p = strrchr(FileName.c_str(),'.');
  return (strstr(p,".Z") == p || strstr(p,".gz") == p );
}


#include <fstream>




StringList SplitString(const string &String, const char Sep)
{
  StringList Result;
  const char* next;
  const char* start = String.c_str();
  if (*start == Sep) ++start;
  while ((next = strchr(start,Sep)))
    {
      char cp[256];
      strncpy(cp, start , next-start);
      cp[next-start] = '\0';
      Result.push_back(string(cp));
      start = next+1;
    }
  Result.push_back(string(start));
  if (getenv("SPLIT_DEBUG"))
    {
      for (StringCIterator p = Result.begin(); p != Result.end(); ++p) cout << *p << ' '; cout << endl;
    }
  return Result;
}


string StandardPath(const string &Path)
/*  Checks that a path does not contain things improperly handled 
    by RelativeLink or correct what can be corrected. */
{
  string standard_path;
  const char *path = Path.c_str();
  if (!path) return NULL;
  if (strstr(path,"./") == path) path +=2;
  standard_path = path;
  if (path[0] == '.' && path[1] == '\0') /* got . */
    {
      //char value[512];
      //getcwd(value, 512);
      standard_path = getenv("PWD");;
    }
  else if (strstr(path,".."))
    {
      string complete_path = FullFileName(path);
      StringList dirs = SplitString(complete_path,'/');
      for (StringIterator p = dirs.begin(); p != dirs.end() ; )
	{
	  if (p->length() == 0) { p = dirs.erase(p); if (p!=dirs.begin()) --p; continue;} 
	  StringIterator next = p;
	  ++next;
	  if (next == dirs.end()) break;
	  if (*p != ".." && *next == "..") // found a/..
	    {
	      p = dirs.erase(p); // erase a
	      p = dirs.erase(p);   // erase ..
	      if (p != dirs.end()) --p; // in case you had b/a/../.., restart at b.
	    }
	  else ++p;
	}
      standard_path = "";
      for (StringCIterator p = dirs.begin(); p != dirs.end() ; ++p)
	{
    
	  standard_path += "/"+(*p);
	}
    }
  return standard_path;
}


/***************************************************************************/
char* StringDuplicate (
		       const char* This
		       )
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  char* string;
  int   length;
  /*.........................................................................*/
  if(This==NULL)   return NULL;
  length           = strlen(This);
  string           = (char*)malloc((length+1)*sizeof(char));
  if(string==NULL) return NULL;
  string[length]   = '\0';
  return           strcpy(string,This);
}


/***************************************************************************/
int StringMatchPattern (
			const char* This   
			,const char* a_pattern 
			)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  int count;
  int      wcount;
  int      lpattern,lstring;
  char*    string;
  char*    token;
  const char* token2;
  int      ok,match;
  int      wnumber;
#define MAX_WORDS 100
  char*    words[MAX_WORDS]; 
  /*.........................................................................*/
  if ( (a_pattern==NULL) && (This==NULL) ) return 1;
  if ( (a_pattern==NULL) && (This!=NULL) ) return 0;
  if ( (a_pattern!=NULL) && (This==NULL) ) return 0;
  lpattern  = strlen(a_pattern);
  lstring   = strlen(This);
  if ((lpattern==0)&&(lstring==0)) return 1;
  if ((lpattern==0)&&(lstring!=0)) return 0;
  if ((lpattern!=0)&&(lstring==0)) return 0;
  /* pattern is * */
  if(strcmp(a_pattern,"*")==0) return 1;
  wcount = 0;
  for(count=0;count<lpattern;count++) {if(a_pattern[count]=='*') wcount++;}
  /* no wildcard */
  if(wcount==0)
    {
      return (strcmp(a_pattern,This)==0 ? 1 : 0 );
    }

  /* complex pattern */
  token    = string = StringDuplicate(a_pattern);
  wnumber  = 0;
  while(1)
    { char* pos;
      pos   = strstr (token,"*");
      if(pos!=NULL)
        {
          *pos = '\0';
          if(*token!='\0') 
            {
              if(wnumber>=MAX_WORDS) {fprintf(stderr, "StringMatchPattern : overflow of MAX_WORDS\n");} 
              else                   {words[wnumber] = token;wnumber++;}
            }
	  token = pos + 1;
        }
      else /*last word*/
        {
          if(*token!='\0') 
            {
              if(wnumber>=MAX_WORDS) {fprintf(stderr, "StringMatchPattern : overflow of MAX_WORDS\n");} 
              else                   {words[wnumber] = token;wnumber++;}
            }
          break;
        }
    }
  /* check that at least one word is not empty */
  ok = 0;
  for(count=0;count<wnumber;count++)
    { 
      if(*(words[count])!='\0') {ok = 1;break;}
    }
  if(ok==0) {free(string);return 1;} /* only wildcards */

  /* loop on words */
  match    = 1;
  token2    = This;
  for(count=0;count<wnumber;count++)
    { int   lword;
      lword = strlen(words[count]); 
      if(lword>0) 
        { 
	  char* pos;
	  if(count==0)
	    {
	      if(a_pattern[0]!='*') /*Begin of pattern (words[0]) and This must match.*/
		{
		  if(strncmp(token2,words[count],lword)!=0) 
		    {
		      match = 0; /*Different.*/
		      break;      
		    }
		  token2 = token2 + lword;
		  continue;
		}
	    }
	  pos                = const_cast<char *>(strstr (token2,words[count]));
	  if(pos==NULL)      {match=0;break;}
	  if((count==(wnumber-1)) && (a_pattern[lpattern-1]!='*') ) /*Last word.*/
	    {
	      if(strcmp(This+lstring-lword,words[count])!=0) match = 0; /*Compare last word and end of This.*/
	      break;
	    }
	  else
	    {
	      token2 = pos + lword;
	    }
        }
    }
  free(string);
  /*printf ("debug:match:%d:%s|%s|\n",match,This,a_pattern);*/
  return        match;
}
/***************************************************************************/
string StringToUpper(const string Source)
//Converts a string to uppercase
{
  int n = Source.length();
  char* dest = new char[n+1];
  for(int i=0; i <= n  ;i++) dest[i] = toupper(Source[i]); 
  return string(dest);
}

string StringToLower(const string Source)
//Converts a string to lowercase
{
  int n = Source.length();
  char* dest = new char[n+1];
  for(int i=0; i <= n  ;i++) dest[i] = tolower(Source[i]); 
  return string(dest);
}

string AddPrefix(const string &Prefix, const string &FileName)
{
  return DirName(FileName)+"/"+Prefix+BaseName(FileName);
}


// fix by Andrew Becker
int RemovePattern(string &Source, const string &aPattern)
{
  // cout << " Removing '" << aPattern << "' from '"<< Source << "'\n";
  int iter = 0;
  size_t pos = Source.find(aPattern);
  if (pos == string::npos) {
    return -1;
  }
  while (pos != (string::size_type) Source.npos)
    {
      Source.erase(pos, aPattern.length());
      iter++;
      pos = Source.find(aPattern);
    }
  return iter;
}

void Replace(string &string_to_be_modified,const string &element_to_replace , const string &element_to_put) {
  string::size_type start = 0;
  string::size_type pos = string_to_be_modified.find(element_to_replace,start);
  while(pos != (string::size_type) string_to_be_modified.npos) {
    string_to_be_modified.replace(pos,element_to_replace.size(),element_to_put,0,element_to_put.size());
    start = pos+element_to_put.size();
    pos = string_to_be_modified.find(element_to_replace,start);
  }  
} 


void DecomposeString(vector<string> &SubStrings, const string &Source,
                     const char *tokens)
{
  string source_copy = Source;
  
  // replace all '\t' by token 
  Replace(source_copy,"\t",tokens);

  string::size_type start = 0;
  string::size_type tokpos = source_copy.find_first_of(tokens,start);
 
  while (tokpos != (string::size_type) source_copy.npos)
    {
      string sub = source_copy.substr(start,tokpos-start);
      RemovePattern(sub," ");
      //RemovePattern(sub,"\t");
      if (sub.length() > 0) SubStrings.push_back(sub);
      start = tokpos+1;
      tokpos = source_copy.find_first_of(tokens,start);
    }
  // get the last element
  string sublast = source_copy.substr(start, source_copy.length());
  RemovePattern(sublast," ");
  //RemovePattern(sublast,"\t");
  if (sublast.length() > 0) SubStrings.push_back(sublast);
                                                                                                               
  // for (size_t i=0; i<SubStrings.size(); ++i)
  //   cout << i << "='"<< SubStrings[i] << "'"<<endl;
}
                                                                                                               
std::string SubstituteExtension(std::string Original, std::string NewExtension)
{
  size_t dotpos = Original.rfind('.');
  if (dotpos > Original.size()) return Original;
  return Original.substr(0,dotpos)+NewExtension;
}

bool CaseInsensitiveCompare(const std::string &s1, const std::string &s2)
{
  int c1, c2;
  const size_t len1(s1.length()), len2(s2.length()), min_len(std::min(len1, len2));
  if (len1 != len2) return false;
  for(size_t i=0; i<min_len; ++i)
    {
      c1=std::tolower(s1[i]);
      c2=std::tolower(s2[i]);
      if(c1!=c2)
	{
	  return false;
	}
    }
  return true;
}
