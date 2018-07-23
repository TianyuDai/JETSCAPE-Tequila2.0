#include "PythiaHadronReader.h"

namespace Jetscape {

template<class T>
PythiaHadronReader<T>::PythiaHadronReader()
{
  VERBOSE(8);
}

template<class T>
PythiaHadronReader<T>::~PythiaHadronReader()
{
  VERBOSE(8);
}

template<class T>
void PythiaHadronReader<T>::Clear()
{
  hadrons.clear();
}

template<class T>
void PythiaHadronReader<T>::AddHadron(string s)
{
	string token;
        strT.set(s);

	vector<string> vS;
        while (!strT.done())
        {
          token = strT.next();
          if(token.compare("H") != 0)
            vS.push_back(token);
        }
        double x[4] = {0.};
        hadrons.push_back(make_shared<Hadron>(stoi(vS[1]),stoi(vS[2]),stoi(vS[3]),stod(vS[4]),stod(vS[5]),stod(vS[6]),stod(vS[7]),x));

}

template<class T>
int PythiaHadronReader<T>::Next(int i)
{
  Clear();

  string line;
  string token;  
  bool flag = false; 
  int position = 0; 
  while (getline(inFile,line))
  {
      strT.set(line); 
      if (strT.isCommentEntry()) continue; 
      if (strT.isEventEntry())
      {
	if (flag) return position; 
	else
	{
          token = strT.next(); 
	  if (i == stoi(token)) flag = true;
	  continue; 
	}
      }
      if (strT.isHadronEntry() && flag) {AddHadron(line); position = inFile.tellg(); }
  }
  return position; 
}

template<class T>
vector<fjcore::PseudoJet>  PythiaHadronReader<T>::GetHadronsForFastJet(){
  vector<fjcore::PseudoJet> forFJ;  

  for ( auto& p : hadrons) {
    forFJ.push_back( p->GetPseudoJet());
  }

  return forFJ;
} 


template<class T>
void PythiaHadronReader<T>::Init()
{
  VERBOSE(8)<<"Open Input File = "<<file_name_in;
  INFO<<"Open Input File = "<<file_name_in;
  
  inFile.open(file_name_in.c_str());
  
  if (!inFile.good())
    { WARN<<"Corrupt input file!"; exit(-1);}
  else
    INFO<<"File opened";

}

template class PythiaHadronReader<ifstream>;

#ifdef USE_GZIP
template class PythiaHadronReader<igzstream>;
#endif

} // end namespace Jetscape
