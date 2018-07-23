#ifndef JETSCAPEPARTONREADER_H
#define JETSCAPEPARTONREADER_H

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
    class JetScapePartonReader
    {

    public:

      JetScapePartonReader();
      JetScapePartonReader(string m_file_name_in) {file_name_in =  m_file_name_in; Init();}
      virtual ~JetScapePartonReader();

      void Close() {inFile.close(); }
      void Clear();
  
      void Next();
      bool Finished() {return inFile.eof();}
      
      int GetCurrentEvent() {return currentEvent-1;}
  
      //shared_ptr<PartonShower> GetPartonShower() {return pShower;}
      //vector<shared_ptr<Parton>> GetPartonList() {return plist;}
      vector<shared_ptr<Parton>> GetPartons() {return partons; }

      vector<fjcore::PseudoJet>  GetPartonsForFastJet();
  
    private:

      StringTokenizer strT;
   
      void Init();
      void AddParton(string s); 
      int currentEvent;
      string file_name_in;
      T inFile;

      //vector<shared_ptr<Parton> > plist;
      vector<shared_ptr<Parton> > partons;
    };

typedef JetScapePartonReader<ifstream> JetScapeReaderAscii;
#ifdef USE_GZIP
typedef JetScapePartonReader<igzstream> JetScapeReaderAsciiGZ;
#endif
  
} // end namespace Jetscape

// ---------------------

#endif

// ---------------------
