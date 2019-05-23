#ifndef _slimmedhit_h_
#define _slimmedhit_h_

#include "TObject.h"

#include<iostream>
#include <vector>

/**
   @short use to strip down the information of a Hit
 */
class SlimmedHit : public TObject
{
 public:
 SlimmedHit() : layerId_(0), assocROIidx_(-1), signalRegion_(-1), en_(0), en_mip_(0), t_(0), cellSize_(0) {}
 SlimmedHit(int layerId, int assocROIidx, int signalRegion, float en, float en_mip, float t, float cellSize) :
  layerId_(layerId), assocROIidx_(assocROIidx), signalRegion_(signalRegion), en_(en), en_mip_(en_mip), t_(t), cellSize_(cellSize) 
  {
  }
  SlimmedHit(const SlimmedHit &other)
    {
      layerId_      = other.layerId_;
      assocROIidx_  = other.assocROIidx_;
      signalRegion_ = other.signalRegion_;      
      en_           = other.en_;
      en_mip_       = other.en_mip_;
      t_            = other.t_;
      cellSize_     = other.cellSize_;
    }

  unsigned int layerId()   { return layerId_; }
  int associatedROI()      { return assocROIidx_; }
  int signalRegion()       { return signalRegion_; }
  float en(bool mip=false) { return mip?en_mip_:en_; }
  float t()                { return t_; } 
  float cellSize()         { return cellSize_; }

  virtual ~SlimmedHit() { }
  
  unsigned int layerId_;
  int   assocROIidx_,signalRegion_;
  float en_,en_mip_,t_;
  float cellSize_;

  ClassDef(SlimmedHit,1)
};

typedef std::vector<SlimmedHit> SlimmedHitCollection;

#endif
