#ifndef DATAFILE__H
#define DATAFILE__H

#include <cstdio>
#include <cstring>
#include <string>

#include "globalval.h"

class DataFile
{
private :
  FILE *f;
  int ncol;
  int ndat;
  std::string filename;
  GlobalVal glob;

  
  public :
    
  DataFile(const std::string &FileName);
  bool IsValid() const {return (f!=NULL);};
  //! number of numerical items per line (assumed to be the same for all lines)
  int NCol() const ;
  //! number of lines of numerical data
  int NDat() const ;

  //!
  const GlobalVal& Globals() const { return glob;}

  //!
  char *next_line(char *line, const int line_length);


  bool next_numerical_data(char* line, const int line_length);
  //! read colums xv=first column, fv: column icol. xv and fv to be new-ed and deleted by caller. return the number of pairs read.
  int ReadCols(double *xv, double *fv, const int ndat, const int Icol=2);

  void Rewind() { if (f) rewind(f);}

  ~DataFile() { if (f) fclose(f);};
};


#endif
