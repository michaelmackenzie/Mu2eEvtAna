//
// Mu2eEvtAna systematics definitions
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_SYSTEMATICS_HH
#define MU2EEVTANA_SYSTEMATICS_HH

//c++ includes
#include <vector>
#include <iostream>

//ROOT includes
#include "TString.h"

//local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"

namespace Mu2eEvtAna {
  class Systematics {
  public:

    // Systematic source (shared between up/down variations)
    enum Type_t {
      kUndefined = -1,
      kNominal   =  0,
      kScale        ,
      kRMCPower
    };

    // Systematic index (histogram index, e.g. obs_<index>)
    enum Sys_t {
      kNominalSys    =  0,
      kScaleUp       =  1,
      kScaleDown     =  2,
      kRMCPowerUp    = 21,
      kRMCPowerDown  = 22
    };

    //-----------------------------------------------------------------------------
    //  data members
    //-----------------------------------------------------------------------------
  public:

    // Define a systematic input
    struct Data_t {
      int     num_;
      Type_t  type_;
      TString name_;
      bool    up_;
      Data_t(int num = -1, Type_t type = kUndefined, TString name = "", bool up = true) : num_(num), type_(type), name_(name), up_(up) {}
    };

    Data_t           data_[kMaxSystematics]; // indexed by systematic number
    std::vector<int> defined_              ; // defined systematic numbers, for looping

    //-----------------------------------------------------------------------------
    //  functions
    //-----------------------------------------------------------------------------
  public:
    Systematics() {
      //initialize the defined systematics information
      for(int isys = 0; isys < kMaxSystematics; ++isys) {
        data_[isys] = getData(isys);
        if(data_[isys].type_ != kUndefined) defined_.push_back(isys); //only loop over defined systematics
      }
    }
    ~Systematics() {}

    const std::vector<int>& Defined() const { return defined_; }

    bool IsDefined(const int isys) const { return isys >= 0 && isys < kMaxSystematics && data_[isys].type_ != kUndefined; }

    Type_t GetType(const int isys) const {
      return (isys >= 0 && isys < kMaxSystematics) ? data_[isys].type_ : kUndefined;
    }

    TString GetName(const int isys) const {
      return (isys >= 0 && isys < kMaxSystematics) ? data_[isys].name_ : TString("");
    }

    int GetNum(const Type_t type, const bool up = true) const {
      for(int isys : defined_) if(data_[isys].type_ == type && data_[isys].up_ == up) return isys;
      std::cout << "Systematics::" << __func__ << ": Undefined systematic type " << type << " (up = " << up << ")" << std::endl;
      return -1;
    }

    int GetNum(const TString name) const {
      for(int isys : defined_) if(data_[isys].up_ && data_[isys].name_ == name) return isys; //map to the up value
      std::cout << "Systematics::" << __func__ << ": Undefined systematic " << name.Data() << std::endl;
      return -1;
    }

    bool IsUp(const int isys) const {
      return (isys >= 0 && isys < kMaxSystematics) ? data_[isys].up_ && data_[isys].type_ != kUndefined : false;
    }
  private:
    static Data_t getData(const int isys) {
      switch(isys) {
      case kNominalSys  : return Data_t(isys, kNominal , "Nominal" , true );
      case kScaleUp     : return Data_t(isys, kScale   , "Scale"   , true );
      case kScaleDown   : return Data_t(isys, kScale   , "Scale"   , false);
      case kRMCPowerUp  : return Data_t(isys, kRMCPower, "RMCPower", true );
      case kRMCPowerDown: return Data_t(isys, kRMCPower, "RMCPower", false);
      default: break;
      }

      //return the default result if no systematic is defined
      return Data_t();
    }
  };
}
#endif
