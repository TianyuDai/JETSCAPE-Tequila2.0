/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>
using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

HardProcess::HardProcess()
{
  VERBOSE(8);
  SetId("HardProcess");
}

HardProcess::~HardProcess()
{
  VERBOSE(8);
  hp_list.clear();
  hd_list.clear();
  disconnect_all();
}

void HardProcess::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize HardProcess : "<<GetId()<< " ...";
 
  fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" );

  if (!fd) {
      WARN << "Not a valid JetScape XML Hard section file or no XML file loaded!";
      exit(-1);
  }
  
  VERBOSE(8);

  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
      WARN << "No initial state module, try: auto trento = make_shared<TrentoInitial>(); jetscape->Add(trento);";
  }
  
  InitTask();

  JetScapeTask::InitTasks();
  plist.push_back(1); 
}

void HardProcess::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();
  
  JetScapeTask::ExecuteTasks();
  plist.push_back(1); 
}

void HardProcess::Clear()
{
  JSDEBUG<<"Clear Hard Process : "<<GetId()<< " ...";

  hp_list.clear();
  hd_list.clear();
  VERBOSE(8)<<hp_list.size();
  //plist.push_back(1); 
}

void HardProcess::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);

  auto f = w.lock();
  if ( f ){
    VERBOSE(8)<<f->GetOutputFileName();
    
    // Weight, xsec, etc
 
    // // Can explicitly write our own header information, though the writer should handle this.
    // std::ostringstream oss;
    // oss.str(""); oss << GetId() << " sigmaGen  = " << GetSigmaGen();  
    // f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " sigmaErr  = " << GetSigmaErr();
    // f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " weight  = " << GetEventWeight();
    // f->WriteComment ( oss.str() );
    
    // Hard partons
    f->WriteComment("HardProcess Parton List: "+GetId());  
    // int i = 0; 
    for ( auto hp : hp_list )
    {
       // f->WriteWhiteSpace("["+to_string(i)+"] P");
       f->Write( hp );
       // ++i; 
    }
    //f->WriteComment("Final State Parton List"); 
  }
}

void HardProcess::CollectHeader( weak_ptr<JetScapeWriter> w ){
  auto f = w.lock();
  if ( f ){
    auto& header = f->GetHeader();
    header.SetSigmaGen( GetSigmaGen() );
    header.SetSigmaErr( GetSigmaErr() );
    header.SetEventWeight( GetEventWeight() );
  }
}

    
} // end namespace Jetscape
