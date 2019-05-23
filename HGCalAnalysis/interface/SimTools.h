#ifndef _simtools_h_
#define _simtools_h_

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

math::XYZVectorD getInteractionPosition(const std::vector<SimTrack> *SimTk, 
                                        const std::vector<SimVertex> *SimVtx, 
                                        int barcode);

#endif
