// $Id: datacards.cc,v 1.8 2011-12-08 11:30:41 betoule Exp $
//
// Datacards, acquisition EROS II
//
//
// Eric Aubourg, Decembre 95
//
// DAPNIA/SPP (Saclay) / CEA

//#include "defs.h"
#include "datacards.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include <cstring> // for strcmp 

using namespace std;

//++
// Class	DataCards
// Lib		Outils++ 
// include	datacards.h
//
//   	Cette classe permet la gestion des parametres d'un programme a partir 
//   	de mot-cle (lecture d'un fichier par exemple)
//--

//++
// Titre	Constructeurs
//--
//++
//
// DataCards()
// DataCards(string const& fn)
//    Createur avec lecture des parametres ds le fichier de nom "fn"
//--

DataCards::DataCards()
{
  fileName = "Nofile";
}

DataCards::DataCards(string const& fn)
{
  fileName = fn;
  ReadFile(fn);
}

//++
// Titre	Methodes
//--
//++
// AddProcF(ProcCard f, string const& mtch="*")
//	Ajoute une fonction de traitement a la liste pour les mots cle
//      compatibles avec la chaine "mtch" ("mtch" peut contenir "*" en debut 
//      fin de mot) 
//
// Clear()
//	Supprime les cartes existantes
//
// ReadFile(string const& fn) 
//      Lit le contenu du fichiers "fn" et ajoute les cartes a la liste
//
// AppendCard(string const& line)
//	Rajoute la carte "line" a la liste
//--

void
DataCards::AddProcF(ProcCard f, string const& mtch)
{
  CrdPF mpf;
  if (f == NULL)  return;
  mpf.pf = f;
  if (mtch.length() <= 0)  mpf.patt = "*";
  else  mpf.patt = mtch;
  cpfs.push_back(mpf);

  // On applique cette nouvelle fonction aux cartes existantes
  CardList::iterator  ic;
  for(ic = cards.begin(); ic != cards.end(); ic ++)
    {
      vector<string>::iterator ik;
      string tks;
      for(ik = (*ic).tokens.begin(); ik != (*ic).tokens.end(); ik++)
	tks = tks + " " + (*ik);
      ApplyPF(mpf, (*ic).kw, tks);
    }
}

void
DataCards::Clear()
{
  cards.erase(cards.begin(), cards.end());
}

void
DataCards::ReadFile(string const& fn)
{
  string file = fn;
  DoReadFile(file);
}

void
DataCards::AppendCard(string const& crd)
{
  Card c;
  size_t p = 1;
  size_t q = crd.find_first_of(" \t");
  size_t l = crd.length();

  string toks;
  if (l < 2)  return;
  if (crd[0] != '@')  return;

  if (q < l)
    {  c.kw = crd.substr(p,q-p);  toks = crd.substr(q, l-q); }
  else { c.kw = crd.substr(p,l-p);  toks = ""; }
  //  On applique les ProcFunc's
  ApplyPFL(c.kw, toks);
  while (q < l) 
    {
      p = crd.find_first_not_of(" \t",q+1); // au debut d'un token
      if (p>=l) break;
      q = crd.find_first_of(" \t",p); // la fin du token;
      string token = crd.substr(p,q-p);
      c.tokens.push_back(token);
    }
  // On supprime la carte de la liste, si elle existe deja ...
  RemoveCard(c.kw);
  cards.push_back(c);
}

void
DataCards::DoReadFile(string const& fn)
{
  char line_buff[1024];
  FILE *fip;
  
  if ( (fip = fopen(fn.c_str(),"r")) == NULL ) 
    {
      cerr << " DataCards::DoReadFile : cannot open " << fn << endl;
      return;
    }
  
  while (fgets(line_buff,1023,fip) != NULL)
    {
      line_buff[strlen(line_buff)-1] = '\0';   /*  LF/CR de la fin */
      string line(line_buff);
      AppendCard(line);
    }
  fclose(fip);
}

int
DataCards::ApplyPF(CrdPF & cpf, string const& key, string const& toks)
{
  int l;
  size_t lk;
  int rc = 0;
  // On verifie si le "pattern" correspond
  bool mtch = false;
  l = cpf.patt.length(); 
  if (cpf.patt == "*")  mtch = true;
  else if (cpf.patt[0] == '*')   
    {
      lk = key.length();
      if (cpf.patt[l-1] != '*')
	{
	  if (strcmp(key.c_str()+(lk-l+1), cpf.patt.c_str()+1) == 0)   mtch = true;
	}
      else if  (key.find(cpf.patt.substr(1,l-2)) < lk)  mtch = true;
    }
  else if (cpf.patt[l-1] == '*')
    {
      if ( strncmp(key.c_str(), cpf.patt.c_str(),l-1) == 0)  mtch = true;
    }
  else if (key == cpf.patt)  mtch = true;

  // Si oui, on appelle la fonction correspondante
  if (mtch)  rc = cpf.pf(key, toks); 

  return(rc);
}


int
DataCards::ApplyPFL(string const& key, string const& toks)
{
  int rc = 0;
  CrdPFList::iterator icf;
  for(icf = cpfs.begin(); icf != cpfs.end(); icf++)
    rc += ApplyPF((*icf), key, toks);
  return(rc);
}

void  
DataCards::RemoveCard(string const& key)
{
  CardList::iterator i;
  for(i=cards.begin(); i != cards.end(); i++)
    if ((*i).kw == key) { cards.erase(i);  break;  }
}

const DataCards::Card *
DataCards::FindKey(string const& key, const bool ThrowIfAbsent) const
{
  /*
    CardList::iterator i = find_if(cards.begin(), cards.end(), bind2nd(KeyEq(),key));
    if (i == cards.end() ) return NULL;
  */
  CardList::const_iterator i;
  for(i=cards.begin(); i != cards.end(); i++)
    if ((*i).kw == key) return &*i;
  if (ThrowIfAbsent)
    {
      cerr << " cannot find key " << key << " in " << fileName << endl;
      throw(MissingKeyException(" cannot find key "+key+" in file "+fileName));
    }
  return NULL;
}

//++
// Titre	Acces aux parametres
//--
//++
// int   NbCards()
//	Renvoie le nombre de cartes data
// bool	 HasKey(string const& key) 
//	Indique l'existence d'une carte avec la cle "key"
// int   NbParam(string const& key)
//	Indique le nombre de parametre (separes par des espaces) pour la cle "key"
// string  SParam(string const& key, int num = 0, string def="")
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) sous forme de
//	chaine de caracteres ("string")
// long    IParam(string const& key, int numero = 0, long def = 0)
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) convertie 
//	en entier ("long")
// double  DParam(string const& key, int numero = 0, double def = 0)
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) convertie 
//	en flottant ("double")
//--


bool
DataCards::HasKey(string const& key) const
{
  return FindKey(key,false) != NULL;
}

int 
DataCards::NbCards() const
{
  return(cards.size());
}

int 
DataCards::NbParam(string const& key) const
{
  const DataCards::Card * p = FindKey(key, true);
  if (!p) return(-1);
  else return(p->tokens.size());
}

string
DataCards::SParam(string const& key, int numero, string def) const
{
  const DataCards::Card * p = FindKey(key, true);
  if (!p) return def;
  if ( (numero < 0) || ((size_t)(numero) >= p->tokens.size()) )  return def;
  return p->tokens[numero];
}

long
DataCards::IParam(string const& key, int numero, long def) const
{
  string p = SParam(key, numero, "");
  if (p == "") return def;
  long i;
  //istrstream(p.c_str(), p.length()) >> i;
  sscanf(p.c_str(),"%ld",&i);
  return i;
}

double
DataCards::DParam(string const& key, int numero, double def) const
{
  string p = SParam(key, numero, "");
  if (p == "") return def;
  double i;
  //istrstream(p.c_str(), p.length()) >> i;
  sscanf(p.c_str(),"%lg",&i);
  return i;
}


vector<string>
DataCards::AllParams(const string &Key) const
{
  const DataCards::Card * p = FindKey(Key, true);
  if (!p) return vector<string>();
  return p->tokens;
}

   
ostream& operator << (ostream& s, const DataCards &c)
{
  for (DataCards::CardList::const_iterator i = c.cards.begin(); i != c.cards.end(); i++) {
    s << setw(10) << (*i).kw << " ";
    for (vector<string>::const_iterator j = (*i).tokens.begin(); j != (*i).tokens.end(); j++)
      s << (*j) << " ";
    s << endl;
  }
  return s;
}



