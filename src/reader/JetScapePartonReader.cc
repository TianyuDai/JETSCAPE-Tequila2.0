#include "JetScapePartonReader.h"

namespace Jetscape {

template<class T>
JetScapePartonReader<T>::JetScapePartonReader()
{
  VERBOSE(8);
  currentEvent=-1;
}

template<class T>
JetScapePartonReader<T>::~JetScapePartonReader()
{
  VERBOSE(8);
}

template<class T>
void JetScapePartonReader<T>::Clear()
{
  partons.clear();
}

template<class T>
void JetScapePartonReader<T>::AddParton(string s)
{
	string token;
        strT.set(s);

	vector<string> vS;
        while (!strT.done())
        {
          token = strT.next();
          if(token.compare("P") != 0)
            vS.push_back(token);
        }
        //double x[4] = {stod(vS[8]), stod(vS[9]), stod(vS[10]), stod(vS[11])};
        double x[4] = {0}; 
        partons.push_back(make_shared<Parton>(stoi(vS[1]),stoi(vS[2]),stoi(vS[3]),stod(vS[4]),stod(vS[5]),stod(vS[6]),stod(vS[7]),x));
}

template<class T>
void JetScapePartonReader<T>::Next()
{
  if (currentEvent>0)
    Clear();

  string line;
  string token;  
  // INFO<<"Current Event = "<<currentEvent;

  //pShowers.push_back(make_shared<PartonShower>());
  //pShower=pShowers[0];
  //currentShower=1;
  
  while (getline(inFile,line))
  {
      strT.set(line); 
      if (strT.isCommentEntry()) continue; 
      if (strT.isEventEntry()) {
	int newEvent=stoi(strT.next());		      
	if (currentEvent!=newEvent && currentEvent>-1) {
	  currentEvent++;
	  break;
	}	
	currentEvent=newEvent;
	continue;
      }
      if (strT.isPartonEntry()) AddParton(line); 
  }
  if (Finished())
    currentEvent++;
}

template<class T>
vector<fjcore::PseudoJet>  JetScapePartonReader<T>::GetPartonsForFastJet(){
  vector<fjcore::PseudoJet> forFJ;  

  for ( auto& p : partons) {
    forFJ.push_back( p->GetPseudoJet());
  }

  return forFJ;
} 


template<class T>
void JetScapePartonReader<T>::Init()
{
  VERBOSE(8)<<"Open Input File = "<<file_name_in;
  INFO<<"Open Input File = "<<file_name_in;
  
  inFile.open(file_name_in.c_str());
  
  if (!inFile.good())
    { WARN<<"Corrupt input file!"; exit(-1);}
  else
    INFO<<"File opened";
    
  currentEvent=0;
}

template class JetScapePartonReader<ifstream>;

#ifdef USE_GZIP
template class JetScapePartonReader<igzstream>;
#endif

} // end namespace Jetscape
