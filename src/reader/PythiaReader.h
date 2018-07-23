#ifndef PYTHIAREADER_H
#define PYTHIAREADER_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
#include "JetScapeParticles.h"
#include "JetScapeLogger.h"
#include "StringTokenizer.h"
#include "PartonShower.h"
#include <fstream>
#ifdef USE_GZIP
#include "gzstream.h"
#endif

using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;

namespace Jetscape {

  template<class T>
    class PythiaReader
    {

    public:

      PythiaReader();
      PythiaReader(string m_file_name_in) {file_name_in =  m_file_name_in; Init();}
      virtual ~PythiaReader();

      void Close() {inFile.close();}
      void Clear();
      void Reset(int p) {inFile.seekg(p);} 
  
      int Next(int i);
      bool Finished() {return inFile.eof();}

      double GetSigmaGenFromFile(int n); 
  
      //shared_ptr<PartonShower> GetPartonShower() {return pShower;}
      vector<shared_ptr<Parton>> GetPartonList() {return plist;}
      vector<shared_ptr<Parton>> GetPartons() {return partons; }

      vector<fjcore::PseudoJet>  GetPartonsForFastJet();
  
    private:

      StringTokenizer strT;
   
      void Init();
      //void AddNode(string s);
      //void AddEdge(string s);
      //void MakeGraph();
      void AddParton(string s); 
      string file_name_in;
      T inFile;
  
      //int currentEvent;
      //int currentShower;
      vector<shared_ptr<Parton> > plist;
      vector<shared_ptr<Parton> > partons;
    };

typedef PythiaReader<ifstream> JetScapeReaderAscii;
#ifdef USE_GZIP
typedef PythiaReader<igzstream> JetScapeReaderAsciiGZ;
#endif
  
} // end namespace Jetscape

// ---------------------

#endif

// ---------------------
