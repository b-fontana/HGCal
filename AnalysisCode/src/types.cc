#include "UserCode/AnalysisCode/interface/types.h"

std::unordered_map<ParticleType, std::string> mParticleType = {{ParticleType::Photon, "Photon"}, 
							       {ParticleType::Pion, "Pion"}};

std::unordered_map<DetectorRegion, std::string> mDetectorRegion = {{DetectorRegion::Inner, "inner"},
								   {DetectorRegion::Central, "central"},
								   {DetectorRegion::Outer, "outer"}};

std::unordered_map<Method, std::string> mMethod= {{Method::ShowerLeakage, "ed"},
						  {Method::BruteForce, "fineeta"}};
							  
