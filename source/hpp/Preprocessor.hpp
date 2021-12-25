#pragma once
#include<iostream>

/**--------------------------------------------*/

#ifdef FVM_NO_EXCEPTIONS
#define FVM_THROW(msg) return false;
#define FVM_COND_THROW(cond, msg) if(cond) return false;
#else
 /// Throws an std::runtime_error with the given message.
#define FVM_THROW(msg) {std::stringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}

/// Throws an std::runtime_error with the given message, if the given condition evaluates to true.
#define FVM_COND_THROW(cond, msg)  if(cond){std::stringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}
#endif

/**--------------------------------------------*/

#define NEW_MACRO(thisClass)									\
  thisClass* thisClass::New() { STANDARD_NEW_BODY(thisClass); }

#define STANDARD_NEW_BODY(thisClass)							\
  thisClass* result = new thisClass;							\
  return result

/**--------------------------------------------*/

#define SET_MACRO(name,type) \
virtual void Set##name (type _arg) \
{ \
    if (this->name != _arg) \
    { \
        this->name = _arg; \
    } \
}

/**--------------------------------------------*/

#define GetMacro(name,type) \
virtual type Get##name () \
{ \
    return this->name; \
}

/**--------------------------------------------*/

#define FREE_OBJ_MACRO(name) \
if(name!=nullptr){  \
    delete name; \
    name = nullptr; \
}

/**--------------------------------------------*/

#define FREE_ARR_MACRO_1P(name) \
if(name!=nullptr){ \
    delete[] name; \
    name = nullptr; \
}

/**--------------------------------------------*/