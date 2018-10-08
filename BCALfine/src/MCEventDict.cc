// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIhenschecdIgeant4dIBCALfinedIsrcdIMCEventDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/henschec/geant4/BCALfine/include/MCEvent.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_fiberHit(void *p = 0);
   static void *newArray_fiberHit(Long_t size, void *p);
   static void delete_fiberHit(void *p);
   static void deleteArray_fiberHit(void *p);
   static void destruct_fiberHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::fiberHit*)
   {
      ::fiberHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::fiberHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("fiberHit", ::fiberHit::Class_Version(), "invalid", 11,
                  typeid(::fiberHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::fiberHit::Dictionary, isa_proxy, 4,
                  sizeof(::fiberHit) );
      instance.SetNew(&new_fiberHit);
      instance.SetNewArray(&newArray_fiberHit);
      instance.SetDelete(&delete_fiberHit);
      instance.SetDeleteArray(&deleteArray_fiberHit);
      instance.SetDestructor(&destruct_fiberHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::fiberHit*)
   {
      return GenerateInitInstanceLocal((::fiberHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::fiberHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_opticalHit(void *p = 0);
   static void *newArray_opticalHit(Long_t size, void *p);
   static void delete_opticalHit(void *p);
   static void deleteArray_opticalHit(void *p);
   static void destruct_opticalHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::opticalHit*)
   {
      ::opticalHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::opticalHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("opticalHit", ::opticalHit::Class_Version(), "invalid", 53,
                  typeid(::opticalHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::opticalHit::Dictionary, isa_proxy, 4,
                  sizeof(::opticalHit) );
      instance.SetNew(&new_opticalHit);
      instance.SetNewArray(&newArray_opticalHit);
      instance.SetDelete(&delete_opticalHit);
      instance.SetDeleteArray(&deleteArray_opticalHit);
      instance.SetDestructor(&destruct_opticalHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::opticalHit*)
   {
      return GenerateInitInstanceLocal((::opticalHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::opticalHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_EventPrimary(void *p = 0);
   static void *newArray_EventPrimary(Long_t size, void *p);
   static void delete_EventPrimary(void *p);
   static void deleteArray_EventPrimary(void *p);
   static void destruct_EventPrimary(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventPrimary*)
   {
      ::EventPrimary *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventPrimary >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EventPrimary", ::EventPrimary::Class_Version(), "invalid", 91,
                  typeid(::EventPrimary), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventPrimary::Dictionary, isa_proxy, 4,
                  sizeof(::EventPrimary) );
      instance.SetNew(&new_EventPrimary);
      instance.SetNewArray(&newArray_EventPrimary);
      instance.SetDelete(&delete_EventPrimary);
      instance.SetDeleteArray(&deleteArray_EventPrimary);
      instance.SetDestructor(&destruct_EventPrimary);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventPrimary*)
   {
      return GenerateInitInstanceLocal((::EventPrimary*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::EventPrimary*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MCEvent(void *p = 0);
   static void *newArray_MCEvent(Long_t size, void *p);
   static void delete_MCEvent(void *p);
   static void deleteArray_MCEvent(void *p);
   static void destruct_MCEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MCEvent*)
   {
      ::MCEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MCEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MCEvent", ::MCEvent::Class_Version(), "invalid", 131,
                  typeid(::MCEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MCEvent::Dictionary, isa_proxy, 4,
                  sizeof(::MCEvent) );
      instance.SetNew(&new_MCEvent);
      instance.SetNewArray(&newArray_MCEvent);
      instance.SetDelete(&delete_MCEvent);
      instance.SetDeleteArray(&deleteArray_MCEvent);
      instance.SetDestructor(&destruct_MCEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MCEvent*)
   {
      return GenerateInitInstanceLocal((::MCEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MCEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MODULE_info(void *p = 0);
   static void *newArray_MODULE_info(Long_t size, void *p);
   static void delete_MODULE_info(void *p);
   static void deleteArray_MODULE_info(void *p);
   static void destruct_MODULE_info(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MODULE_info*)
   {
      ::MODULE_info *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MODULE_info >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MODULE_info", ::MODULE_info::Class_Version(), "invalid", 182,
                  typeid(::MODULE_info), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MODULE_info::Dictionary, isa_proxy, 4,
                  sizeof(::MODULE_info) );
      instance.SetNew(&new_MODULE_info);
      instance.SetNewArray(&newArray_MODULE_info);
      instance.SetDelete(&delete_MODULE_info);
      instance.SetDeleteArray(&deleteArray_MODULE_info);
      instance.SetDestructor(&destruct_MODULE_info);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MODULE_info*)
   {
      return GenerateInitInstanceLocal((::MODULE_info*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MODULE_info*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr fiberHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *fiberHit::Class_Name()
{
   return "fiberHit";
}

//______________________________________________________________________________
const char *fiberHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fiberHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int fiberHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fiberHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *fiberHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fiberHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *fiberHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fiberHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr opticalHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *opticalHit::Class_Name()
{
   return "opticalHit";
}

//______________________________________________________________________________
const char *opticalHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::opticalHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int opticalHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::opticalHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *opticalHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::opticalHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *opticalHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::opticalHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EventPrimary::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventPrimary::Class_Name()
{
   return "EventPrimary";
}

//______________________________________________________________________________
const char *EventPrimary::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventPrimary*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventPrimary::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventPrimary*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventPrimary::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventPrimary*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventPrimary::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventPrimary*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MCEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCEvent::Class_Name()
{
   return "MCEvent";
}

//______________________________________________________________________________
const char *MCEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MCEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MCEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MCEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MCEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MODULE_info::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MODULE_info::Class_Name()
{
   return "MODULE_info";
}

//______________________________________________________________________________
const char *MODULE_info::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MODULE_info*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MODULE_info::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MODULE_info*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MODULE_info::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MODULE_info*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MODULE_info::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MODULE_info*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void fiberHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class fiberHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(fiberHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(fiberHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_fiberHit(void *p) {
      return  p ? new(p) ::fiberHit : new ::fiberHit;
   }
   static void *newArray_fiberHit(Long_t nElements, void *p) {
      return p ? new(p) ::fiberHit[nElements] : new ::fiberHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_fiberHit(void *p) {
      delete ((::fiberHit*)p);
   }
   static void deleteArray_fiberHit(void *p) {
      delete [] ((::fiberHit*)p);
   }
   static void destruct_fiberHit(void *p) {
      typedef ::fiberHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::fiberHit

//______________________________________________________________________________
void opticalHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class opticalHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(opticalHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(opticalHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_opticalHit(void *p) {
      return  p ? new(p) ::opticalHit : new ::opticalHit;
   }
   static void *newArray_opticalHit(Long_t nElements, void *p) {
      return p ? new(p) ::opticalHit[nElements] : new ::opticalHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_opticalHit(void *p) {
      delete ((::opticalHit*)p);
   }
   static void deleteArray_opticalHit(void *p) {
      delete [] ((::opticalHit*)p);
   }
   static void destruct_opticalHit(void *p) {
      typedef ::opticalHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::opticalHit

//______________________________________________________________________________
void EventPrimary::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventPrimary.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EventPrimary::Class(),this);
   } else {
      R__b.WriteClassBuffer(EventPrimary::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventPrimary(void *p) {
      return  p ? new(p) ::EventPrimary : new ::EventPrimary;
   }
   static void *newArray_EventPrimary(Long_t nElements, void *p) {
      return p ? new(p) ::EventPrimary[nElements] : new ::EventPrimary[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventPrimary(void *p) {
      delete ((::EventPrimary*)p);
   }
   static void deleteArray_EventPrimary(void *p) {
      delete [] ((::EventPrimary*)p);
   }
   static void destruct_EventPrimary(void *p) {
      typedef ::EventPrimary current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EventPrimary

//______________________________________________________________________________
void MCEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class MCEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MCEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(MCEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MCEvent(void *p) {
      return  p ? new(p) ::MCEvent : new ::MCEvent;
   }
   static void *newArray_MCEvent(Long_t nElements, void *p) {
      return p ? new(p) ::MCEvent[nElements] : new ::MCEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_MCEvent(void *p) {
      delete ((::MCEvent*)p);
   }
   static void deleteArray_MCEvent(void *p) {
      delete [] ((::MCEvent*)p);
   }
   static void destruct_MCEvent(void *p) {
      typedef ::MCEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MCEvent

//______________________________________________________________________________
void MODULE_info::Streamer(TBuffer &R__b)
{
   // Stream an object of class MODULE_info.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MODULE_info::Class(),this);
   } else {
      R__b.WriteClassBuffer(MODULE_info::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MODULE_info(void *p) {
      return  p ? new(p) ::MODULE_info : new ::MODULE_info;
   }
   static void *newArray_MODULE_info(Long_t nElements, void *p) {
      return p ? new(p) ::MODULE_info[nElements] : new ::MODULE_info[nElements];
   }
   // Wrapper around operator delete
   static void delete_MODULE_info(void *p) {
      delete ((::MODULE_info*)p);
   }
   static void deleteArray_MODULE_info(void *p) {
      delete [] ((::MODULE_info*)p);
   }
   static void destruct_MODULE_info(void *p) {
      typedef ::MODULE_info current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MODULE_info

namespace ROOT {
   static TClass *vectorlEpairlEdoublecOdoublegRsPgR_Dictionary();
   static void vectorlEpairlEdoublecOdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEdoublecOdoublegRsPgR(void *p = 0);
   static void *newArray_vectorlEpairlEdoublecOdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEdoublecOdoublegRsPgR(void *p);
   static void deleteArray_vectorlEpairlEdoublecOdoublegRsPgR(void *p);
   static void destruct_vectorlEpairlEdoublecOdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<double,double> >*)
   {
      vector<pair<double,double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<double,double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<double,double> >", -2, "vector", 214,
                  typeid(vector<pair<double,double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEdoublecOdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<double,double> >) );
      instance.SetNew(&new_vectorlEpairlEdoublecOdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEdoublecOdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEdoublecOdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEdoublecOdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEdoublecOdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<double,double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<pair<double,double> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEdoublecOdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<double,double> >*)0x0)->GetClass();
      vectorlEpairlEdoublecOdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEdoublecOdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEdoublecOdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<double,double> > : new vector<pair<double,double> >;
   }
   static void *newArray_vectorlEpairlEdoublecOdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<double,double> >[nElements] : new vector<pair<double,double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEdoublecOdoublegRsPgR(void *p) {
      delete ((vector<pair<double,double> >*)p);
   }
   static void deleteArray_vectorlEpairlEdoublecOdoublegRsPgR(void *p) {
      delete [] ((vector<pair<double,double> >*)p);
   }
   static void destruct_vectorlEpairlEdoublecOdoublegRsPgR(void *p) {
      typedef vector<pair<double,double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<double,double> >

namespace {
  void TriggerDictionaryInitialization_MCEventDict_Impl() {
    static const char* headers[] = {
"/home/henschec/geant4/BCALfine/include/MCEvent.hh",
0
    };
    static const char* includePaths[] = {
"/opt/ROOT/root_6.06.06-install/include/root",
"/home/henschec/geant4/BCALfine_build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MCEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP([Analyze] Fiber Event /*do not remove this comment*/)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/henschec/geant4/BCALfine/include/MCEvent.hh")))  fiberHit;
class __attribute__((annotate(R"ATTRDUMP([Analyze] Fiber Event /*do not remove this comment*/)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/henschec/geant4/BCALfine/include/MCEvent.hh")))  opticalHit;
class __attribute__((annotate(R"ATTRDUMP([Analyze] event primary /*do not remove this comment*/)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/henschec/geant4/BCALfine/include/MCEvent.hh")))  EventPrimary;
class __attribute__((annotate(R"ATTRDUMP([Analyze] MCEvent structure /*do not remove this comment*/)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/henschec/geant4/BCALfine/include/MCEvent.hh")))  MCEvent;
class __attribute__((annotate("$clingAutoload$/home/henschec/geant4/BCALfine/include/MCEvent.hh")))  MODULE_info;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MCEventDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/henschec/geant4/BCALfine/include/MCEvent.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"EventPrimary", payloadCode, "@",
"MCEvent", payloadCode, "@",
"MODULE_info", payloadCode, "@",
"fiberHit", payloadCode, "@",
"opticalHit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MCEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MCEventDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MCEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MCEventDict() {
  TriggerDictionaryInitialization_MCEventDict_Impl();
}
