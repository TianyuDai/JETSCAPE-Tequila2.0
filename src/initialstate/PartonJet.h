#ifndef PARTONJET_H
#define PARTONJET_H
#include "JetClass.h"
#include <vector>

namespace Jetscape{
class PartonJet : public JetScapeModuleBase
{
  public: 
    PartonJet(); 
    virtual ~PartonJet(); 
    virtual void Init(); 
    virtual void Exec(); 
    virtual void Clear(); 
    void AddParton(); 
}; 
}
#endif
