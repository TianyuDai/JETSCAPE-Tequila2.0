#ifndef PYTHIAHADRONREADER_H
#define PYTHIARHADRONEADER_H

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
    class PythiaHadronReader
    {

    public:

      PythiaHadronReader();
      PythiaHadronReader(string m_file_name_in) {file_name_in =  m_file_name_in; Init();}
      virtual ~PythiaHadronReader();

      void Close() {inFile.close();}
      void Clear();
      void Reset(int p) {inFile.seekg(p);} 
  
      int Next(int i);
      bool Finished() {return inFile.eof();}

      vector<shared_ptr<Hadron>> GetHadrons() {return hadrons; }

      vector<fjcore::PseudoJet>  GetHadronsForFastJet();
  
    private:

      StringTokenizer strT;
   
      void Init();
      void AddHadron(string s); 
      string file_name_in;
      T inFile;

      vector<shared_ptr<Hadron> > hadrons;
    };

typedef PythiaHadronReader<ifstream> JetScapeReaderAscii;
#ifdef USE_GZIP
typedef PythiaHadronReader<igzstream> JetScapeReaderAsciiGZ;
#endif
  
} // end namespace Jetscape

// ---------------------

#endif

// ---------------------
