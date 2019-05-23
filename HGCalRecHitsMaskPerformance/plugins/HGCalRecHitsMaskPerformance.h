#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using reco::TrackCollection;

class HGCalRecHitsMaskPerformance : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HGCalRecHitsMaskPerformance(const edm::ParameterSet&);
      ~HGCalRecHitsMaskPerformance();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //variables
      edm::EDGetTokenT<TrackCollection> tracksToken_;  
};
