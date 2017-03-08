#ifndef PATHTOSIMUL__SEEN
#define PATHTOSIMUL__SEEN
#include <string>

/*! returns directory where data are (set by variable SALTPATH, or old PATHMODEL)*/
std::string GetDataPath(); 
std::string DirName(const std::string &FileName);
std::string AssembleFileName(const std::string &DirName, const std::string &FileName, 
			const bool Warn = true);

std::string LocateFile(const std::string &FileName);
std::string LocateFileFromCard(const std::string &Key);

#endif /* PATHTOSIMUL__SEEN*/
