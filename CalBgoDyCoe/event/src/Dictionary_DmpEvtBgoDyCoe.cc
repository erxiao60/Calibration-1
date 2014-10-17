//
// File generated by rootcint at Tue Oct 14 19:12:00 2014

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dOdIincludedIDictionary_DmpEvtBgoDyCoe
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "Dictionary_DmpEvtBgoDyCoe.h"

#include "TCollectionProxyInfo.h"
#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void DmpEvtBgoDyCoe_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_DmpEvtBgoDyCoe(void *p = 0);
   static void *newArray_DmpEvtBgoDyCoe(Long_t size, void *p);
   static void delete_DmpEvtBgoDyCoe(void *p);
   static void deleteArray_DmpEvtBgoDyCoe(void *p);
   static void destruct_DmpEvtBgoDyCoe(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DmpEvtBgoDyCoe*)
   {
      ::DmpEvtBgoDyCoe *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DmpEvtBgoDyCoe >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DmpEvtBgoDyCoe", ::DmpEvtBgoDyCoe::Class_Version(), "./include/DmpEvtBgoDyCoe.h", 13,
                  typeid(::DmpEvtBgoDyCoe), DefineBehavior(ptr, ptr),
                  &::DmpEvtBgoDyCoe::Dictionary, isa_proxy, 4,
                  sizeof(::DmpEvtBgoDyCoe) );
      instance.SetNew(&new_DmpEvtBgoDyCoe);
      instance.SetNewArray(&newArray_DmpEvtBgoDyCoe);
      instance.SetDelete(&delete_DmpEvtBgoDyCoe);
      instance.SetDeleteArray(&deleteArray_DmpEvtBgoDyCoe);
      instance.SetDestructor(&destruct_DmpEvtBgoDyCoe);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DmpEvtBgoDyCoe*)
   {
      return GenerateInitInstanceLocal((::DmpEvtBgoDyCoe*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DmpEvtBgoDyCoe*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *DmpEvtBgoDyCoe::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *DmpEvtBgoDyCoe::Class_Name()
{
   return "DmpEvtBgoDyCoe";
}

//______________________________________________________________________________
const char *DmpEvtBgoDyCoe::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmpEvtBgoDyCoe*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DmpEvtBgoDyCoe::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmpEvtBgoDyCoe*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void DmpEvtBgoDyCoe::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmpEvtBgoDyCoe*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *DmpEvtBgoDyCoe::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmpEvtBgoDyCoe*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void DmpEvtBgoDyCoe::Streamer(TBuffer &R__b)
{
   // Stream an object of class DmpEvtBgoDyCoe.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(DmpEvtBgoDyCoe::Class(),this);
   } else {
      R__b.WriteClassBuffer(DmpEvtBgoDyCoe::Class(),this);
   }
}

//______________________________________________________________________________
void DmpEvtBgoDyCoe::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class DmpEvtBgoDyCoe.
      TClass *R__cl = ::DmpEvtBgoDyCoe::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "UsedFileName", (void*)&UsedFileName);
      R__insp.InspectMember("string", (void*)&UsedFileName, "UsedFileName.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "StartTime", &StartTime);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "StopTime", &StopTime);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "GlobalPMTID", (void*)&GlobalPMTID);
      R__insp.InspectMember("vector<double>", (void*)&GlobalPMTID, "GlobalPMTID.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Inc_Dy5vsDy2", (void*)&Inc_Dy5vsDy2);
      R__insp.InspectMember("vector<double>", (void*)&Inc_Dy5vsDy2, "Inc_Dy5vsDy2.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Inc_Dy8vsDy5", (void*)&Inc_Dy8vsDy5);
      R__insp.InspectMember("vector<double>", (void*)&Inc_Dy8vsDy5, "Inc_Dy8vsDy5.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Slp_Dy5vsDy2", (void*)&Slp_Dy5vsDy2);
      R__insp.InspectMember("vector<double>", (void*)&Slp_Dy5vsDy2, "Slp_Dy5vsDy2.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Slp_Dy8vsDy5", (void*)&Slp_Dy8vsDy5);
      R__insp.InspectMember("vector<double>", (void*)&Slp_Dy8vsDy5, "Slp_Dy8vsDy5.", false);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DmpEvtBgoDyCoe(void *p) {
      return  p ? new(p) ::DmpEvtBgoDyCoe : new ::DmpEvtBgoDyCoe;
   }
   static void *newArray_DmpEvtBgoDyCoe(Long_t nElements, void *p) {
      return p ? new(p) ::DmpEvtBgoDyCoe[nElements] : new ::DmpEvtBgoDyCoe[nElements];
   }
   // Wrapper around operator delete
   static void delete_DmpEvtBgoDyCoe(void *p) {
      delete ((::DmpEvtBgoDyCoe*)p);
   }
   static void deleteArray_DmpEvtBgoDyCoe(void *p) {
      delete [] ((::DmpEvtBgoDyCoe*)p);
   }
   static void destruct_DmpEvtBgoDyCoe(void *p) {
      typedef ::DmpEvtBgoDyCoe current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DmpEvtBgoDyCoe

namespace ROOT {
   void vectorlEdoublegR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEdoublegR_Dictionary();
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector.dll", 0,
                  typeid(vector<double>), DefineBehavior(ptr, ptr),
                  0, &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEdoublegR_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

/********************************************************
* ./include/Dictionary_DmpEvtBgoDyCoe.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableDictionary_DmpEvtBgoDyCoe();

extern "C" void G__set_cpp_environmentDictionary_DmpEvtBgoDyCoe() {
  G__cpp_reset_tagtableDictionary_DmpEvtBgoDyCoe();
}
#include <new>
extern "C" int G__cpp_dllrevDictionary_DmpEvtBgoDyCoe() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmpEvtBgoDyCoe */
static int G__Dictionary_DmpEvtBgoDyCoe_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmpEvtBgoDyCoe* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmpEvtBgoDyCoe[n];
     } else {
       p = new((void*) gvp) DmpEvtBgoDyCoe[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmpEvtBgoDyCoe;
     } else {
       p = new((void*) gvp) DmpEvtBgoDyCoe;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmpEvtBgoDyCoe* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new DmpEvtBgoDyCoe(*(DmpEvtBgoDyCoe*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) DmpEvtBgoDyCoe(*(DmpEvtBgoDyCoe*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmpEvtBgoDyCoe* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new DmpEvtBgoDyCoe(libp->para[0].ref ? *(const DmpEvtBgoDyCoe**) libp->para[0].ref : *(const DmpEvtBgoDyCoe**) (void*) (&G__Mlong(libp->para[0])));
   } else {
     p = new((void*) gvp) DmpEvtBgoDyCoe(libp->para[0].ref ? *(const DmpEvtBgoDyCoe**) libp->para[0].ref : *(const DmpEvtBgoDyCoe**) (void*) (&G__Mlong(libp->para[0])));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmpEvtBgoDyCoe*) G__getstructoffset())->Reset();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) DmpEvtBgoDyCoe::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmpEvtBgoDyCoe::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) DmpEvtBgoDyCoe::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmpEvtBgoDyCoe::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmpEvtBgoDyCoe*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmpEvtBgoDyCoe::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmpEvtBgoDyCoe::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmpEvtBgoDyCoe::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Dictionary_DmpEvtBgoDyCoe_168_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmpEvtBgoDyCoe::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef DmpEvtBgoDyCoe G__TDmpEvtBgoDyCoe;
static int G__Dictionary_DmpEvtBgoDyCoe_168_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (DmpEvtBgoDyCoe*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((DmpEvtBgoDyCoe*) (soff+(sizeof(DmpEvtBgoDyCoe)*i)))->~G__TDmpEvtBgoDyCoe();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (DmpEvtBgoDyCoe*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((DmpEvtBgoDyCoe*) (soff))->~G__TDmpEvtBgoDyCoe();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__Dictionary_DmpEvtBgoDyCoe_168_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmpEvtBgoDyCoe* dest = (DmpEvtBgoDyCoe*) G__getstructoffset();
   *dest = *(DmpEvtBgoDyCoe*) libp->para[0].ref;
   const DmpEvtBgoDyCoe& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmpEvtBgoDyCoe */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDictionary_DmpEvtBgoDyCoe {
 public:
  G__Sizep2memfuncDictionary_DmpEvtBgoDyCoe(): p(&G__Sizep2memfuncDictionary_DmpEvtBgoDyCoe::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDictionary_DmpEvtBgoDyCoe::*p)();
};

size_t G__get_sizep2memfuncDictionary_DmpEvtBgoDyCoe()
{
  G__Sizep2memfuncDictionary_DmpEvtBgoDyCoe a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceDictionary_DmpEvtBgoDyCoe() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe))) {
     DmpEvtBgoDyCoe *G__Lderived;
     G__Lderived=(DmpEvtBgoDyCoe*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe),G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDictionary_DmpEvtBgoDyCoe() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmpEvtBgoDyCoe */
static void G__setup_memvarDmpEvtBgoDyCoe(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe));
   { DmpEvtBgoDyCoe *p; p=(DmpEvtBgoDyCoe*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->UsedFileName)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_string),-1,-1,1,"UsedFileName=",0,"file name of raw data");
   G__memvar_setup((void*)((long)(&p->StartTime)-(long)(p)),105,0,0,-1,-1,-1,1,"StartTime=",0,"the time of the first event used to cal . DyCoe");
   G__memvar_setup((void*)((long)(&p->StopTime)-(long)(p)),105,0,0,-1,-1,-1,1,"StopTime=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->GlobalPMTID)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR),G__defined_typename("vector<double>"),-1,1,"GlobalPMTID=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Inc_Dy5vsDy2)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR),G__defined_typename("vector<double>"),-1,1,"Inc_Dy5vsDy2=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Inc_Dy8vsDy5)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR),G__defined_typename("vector<double>"),-1,1,"Inc_Dy8vsDy5=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Slp_Dy5vsDy2)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR),G__defined_typename("vector<double>"),-1,1,"Slp_Dy5vsDy2=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Slp_Dy8vsDy5)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR),G__defined_typename("vector<double>"),-1,1,"Slp_Dy8vsDy5=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDictionary_DmpEvtBgoDyCoe() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncDmpEvtBgoDyCoe(void) {
   /* DmpEvtBgoDyCoe */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe));
   G__memfunc_setup("DmpEvtBgoDyCoe",1340,G__Dictionary_DmpEvtBgoDyCoe_168_0_1, 105, G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DmpEvtBgoDyCoe",1340,G__Dictionary_DmpEvtBgoDyCoe_168_0_2, 105, G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe), -1, 0, 1, 1, 1, 0, "u 'DmpEvtBgoDyCoe' - 11 - r", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DmpEvtBgoDyCoe",1340,G__Dictionary_DmpEvtBgoDyCoe_168_0_3, 105, G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe), -1, 0, 1, 1, 1, 0, "U 'DmpEvtBgoDyCoe' - 11 - r", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Reset",515,G__Dictionary_DmpEvtBgoDyCoe_168_0_4, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__Dictionary_DmpEvtBgoDyCoe_168_0_5, 85, G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&DmpEvtBgoDyCoe::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__Dictionary_DmpEvtBgoDyCoe_168_0_6, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmpEvtBgoDyCoe::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__Dictionary_DmpEvtBgoDyCoe_168_0_7, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&DmpEvtBgoDyCoe::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__Dictionary_DmpEvtBgoDyCoe_168_0_8, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&DmpEvtBgoDyCoe::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__Dictionary_DmpEvtBgoDyCoe_168_0_12, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__Dictionary_DmpEvtBgoDyCoe_168_0_13, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmpEvtBgoDyCoe::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__Dictionary_DmpEvtBgoDyCoe_168_0_14, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmpEvtBgoDyCoe::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__Dictionary_DmpEvtBgoDyCoe_168_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmpEvtBgoDyCoe::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__Dictionary_DmpEvtBgoDyCoe_168_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmpEvtBgoDyCoe::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~DmpEvtBgoDyCoe", 1466, G__Dictionary_DmpEvtBgoDyCoe_168_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__Dictionary_DmpEvtBgoDyCoe_168_0_18, (int) ('u'), G__get_linked_tagnum(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe), -1, 1, 1, 1, 1, 0, "u 'DmpEvtBgoDyCoe' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDictionary_DmpEvtBgoDyCoe() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalDictionary_DmpEvtBgoDyCoe() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcDictionary_DmpEvtBgoDyCoe() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR = { "vector<double,allocator<double> >" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe = { "DmpEvtBgoDyCoe" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDictionary_DmpEvtBgoDyCoe() {
  G__Dictionary_DmpEvtBgoDyCoeLN_TClass.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_TBuffer.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_TMemberInspector.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_TObject.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_string.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDictionary_DmpEvtBgoDyCoe() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_TClass);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_TObject);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEdoublecOallocatorlEdoublegRsPgR);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_string);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__Dictionary_DmpEvtBgoDyCoeLN_DmpEvtBgoDyCoe),sizeof(DmpEvtBgoDyCoe),-1,292608,(char*)NULL,G__setup_memvarDmpEvtBgoDyCoe,G__setup_memfuncDmpEvtBgoDyCoe);
}
extern "C" void G__cpp_setupDictionary_DmpEvtBgoDyCoe(void) {
  G__check_setup_version(30051515,"G__cpp_setupDictionary_DmpEvtBgoDyCoe()");
  G__set_cpp_environmentDictionary_DmpEvtBgoDyCoe();
  G__cpp_setup_tagtableDictionary_DmpEvtBgoDyCoe();

  G__cpp_setup_inheritanceDictionary_DmpEvtBgoDyCoe();

  G__cpp_setup_typetableDictionary_DmpEvtBgoDyCoe();

  G__cpp_setup_memvarDictionary_DmpEvtBgoDyCoe();

  G__cpp_setup_memfuncDictionary_DmpEvtBgoDyCoe();
  G__cpp_setup_globalDictionary_DmpEvtBgoDyCoe();
  G__cpp_setup_funcDictionary_DmpEvtBgoDyCoe();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDictionary_DmpEvtBgoDyCoe();
  return;
}
class G__cpp_setup_initDictionary_DmpEvtBgoDyCoe {
  public:
    G__cpp_setup_initDictionary_DmpEvtBgoDyCoe() { G__add_setup_func("Dictionary_DmpEvtBgoDyCoe",(G__incsetup)(&G__cpp_setupDictionary_DmpEvtBgoDyCoe)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDictionary_DmpEvtBgoDyCoe() { G__remove_setup_func("Dictionary_DmpEvtBgoDyCoe"); }
};
G__cpp_setup_initDictionary_DmpEvtBgoDyCoe G__cpp_setup_initializerDictionary_DmpEvtBgoDyCoe;
