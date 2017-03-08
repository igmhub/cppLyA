#include <iostream>
#include "datafile.h"
#include <string>

#include <ctype.h> // for isdigit and is alpha */
#include <cstdlib> // for strtod

#include "snfitexception.h"

using namespace std;


bool DataFile::next_numerical_data(char* line, const int line_length)
{
  char *where;
  do 
    where = next_line(line, line_length); 
  while (where &&  !(strpbrk(where,"0123456789.+- ") == where));
  return (where != NULL);
}



char *DataFile::next_line(char *line, const int line_length)
{
  char *p = NULL;
   while ((p = fgets(line, line_length,f)))
    {
      // skip comments
      if (line[0] == '#') continue;
      // lines with '@' are read otherwise 
      if (line[0] == '@') continue;
      // skip blank lines
      p = line;
      while (*p == ' ' || *p == '\t') ++p;
      if (*p == '\n' || *p == '\0') continue;
      break;
    }
  return p;
}
  

DataFile::DataFile(const string &FileName)
{
  f = fopen(FileName.c_str(),"r");
  filename = FileName;
  ncol = 0;
  if (!f)
    {
      cerr << " Cannot open " << FileName << endl;
      throw(SnfitException(" DataFile : cannot open \""+FileName+"\""));
    }

  char line[1024];
  while (fgets(line, 1024,f))
    if (line[0] == '@') glob.ProcessLine(line);
  rewind(f);

  ndat = 0;

  if (next_numerical_data(line,1024) == false) 
    {
      //      cerr << " no numerical data found in " << FileName << endl;
      rewind(f);
      return;
    }  
  // count the number of (numerical) items per lines
  char *p = line;
  char *next = p;
  ncol=0; 
  for (;;) 
    { 
      strtod(p, &next);
      if (next == p) break;
      p = next;
      ncol++;
    }

  // count the number of (numerical) lines
  ndat = 1; // because we alredy read 1 ( to count the number of columns)
  while (next_numerical_data(line,1024)) ndat++;
  rewind(f);

  // Load the "GlobalVals (@KEY Value(s) lines)
}

int DataFile::NCol() const
{
  return ncol;
}


int DataFile::NDat() const
{
  return ndat;
}



static int split_numerical_line(char *line, double *values)
{
  char *p, *next;
  p = line;
  int count=0;
  for (;;)
    {
      values[count] = strtod(p,&next);
      if (next == p) break;
      p = next;
      count++;
    }
  return count;
}




int DataFile::ReadCols(double *xv, double *fv, const int ndat, const int Icol)
{ 
  char line[512];
  int count = 0;
  double values[50];
  if (f) rewind(f); else return 0;
  while((next_numerical_data(line, 512)))
    {
      if (count == ndat)
	{
	  count ++; // to trigger the error message
          break;
	}
      if (split_numerical_line(line,values) < Icol)
	{
	  cerr << " ERROR : not enough columns in " << filename  << " to read column " << Icol << endl;
	  cerr << " ERROR : faulty line follows :"  << endl << line << endl;
	  break;
	}
      xv[count] = values[0];
      fv[count] = values[Icol-1];
      count ++;
    }
  if (count != ndat)
    {
      cerr << " General1DFunction::General1DFunction : Inconsistency between DataFile::NDat() and " << endl 
	   << " ...... just counting the numerical data line in " << filename << " end of file not read ... " << endl;
    }
  return count;
}
