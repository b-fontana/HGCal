#include "UserCode/HGCalAnalysis/interface/SimTools.h"

#include <iostream>

using namespace std;

//
math::XYZVectorD getInteractionPosition(const std::vector<SimTrack> *SimTk, 
                                        const std::vector<SimVertex> *SimVtx, 
                                        int barcode)
{
  math::XYZVectorD x(0,0,0);
  
  //loop over vertices
  for (const SimVertex &simVtx : *SimVtx) 
    {
      //require the parent to be the given barcode
      bool noParent( simVtx.noParent() );
      if(noParent) continue;
      int pIdx( simVtx.parentIndex() );
      if( pIdx!=barcode) continue;      
      x=math::XYZVectorD(simVtx.position());
      break;
    }

  return x;
}
