#include <iostream>
#include <vector>
#include <iterator>
#include "TCanvas.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "TProfile.h"

#include "interface/utils.h"
#include "interface/calibration.h"
#include "interface/software_correction.h"

int_ main() {
  unsigned int ncores = std::thread::hardware_concurrency();

  //variables
  float_ mingenen = 20.;
  vec1d<float_> etareg = {2.7, 2.94};
  vec1d<float_> etareg_shift = VecOps(etareg).shift();
  uint_ n = etareg.size();
  const int_ nreg = 3;
  const int_ nlayers = 28;
  std::string label = "";
  std::string samples = "inner";
  uint_ mask = 3;
  std::string noPUFile = "../summaries/summary_mask3_inner.root";
  std::string outpath = "../pics_inner/mask3";
  vec1d<float_> bckgcuts = {0.5, 0.6, 0.7};

  std::string def1 = 
    "std::vector<float> en = {en_sr1_ROI, en_sr2_ROI, en_sr3_ROI};"
    "return en;";
  std::string def2 = 
    "std::vector<float> noise = {noise_sr1_ROI, noise_sr2_ROI, noise_sr3_ROI};"
    "return noise;";
  std::string def3 = 
    "std::vector<float> en_layer1 = {en_sr1_layer1, en_sr1_layer2, en_sr1_layer3,"
    "en_sr1_layer4, en_sr1_layer5, en_sr1_layer6, en_sr1_layer7, en_sr1_layer8,"
    "en_sr1_layer9, en_sr1_layer10, en_sr1_layer11, en_sr1_layer12, en_sr1_layer13,"
    "en_sr1_layer14, en_sr1_layer15, en_sr1_layer16, en_sr1_layer17, en_sr1_layer18,"
    "en_sr1_layer19, en_sr1_layer20, en_sr1_layer21, en_sr1_layer22, en_sr1_layer23,"
    "en_sr1_layer24, en_sr1_layer25, en_sr1_layer26, en_sr1_layer27, en_sr1_layer28};"
    "std::vector<float> en_layer2 = {en_sr2_layer1, en_sr2_layer2, en_sr2_layer3,"
    "en_sr2_layer4, en_sr2_layer5, en_sr2_layer6, en_sr2_layer7, en_sr2_layer8,"
    "en_sr2_layer9, en_sr2_layer10, en_sr2_layer11, en_sr2_layer12, en_sr2_layer13,"
    "en_sr2_layer14, en_sr2_layer15, en_sr2_layer16, en_sr2_layer17, en_sr2_layer18,"
    "en_sr2_layer19, en_sr2_layer20, en_sr2_layer21, en_sr2_layer22, en_sr2_layer23,"
    "en_sr2_layer24, en_sr2_layer25, en_sr2_layer26, en_sr2_layer27, en_sr2_layer28};"
    "std::vector<float> en_layer3 = {en_sr3_layer1, en_sr3_layer2, en_sr3_layer3,"
    "en_sr3_layer4, en_sr3_layer5, en_sr3_layer6, en_sr3_layer7, en_sr3_layer8,"
    "en_sr3_layer9, en_sr3_layer10, en_sr3_layer11, en_sr3_layer12, en_sr3_layer13,"
    "en_sr3_layer14, en_sr3_layer15, en_sr3_layer16, en_sr3_layer17, en_sr3_layer18,"
    "en_sr3_layer19, en_sr3_layer20, en_sr3_layer21, en_sr3_layer22, en_sr3_layer23,"
    "en_sr3_layer24, en_sr3_layer25, en_sr3_layer26, en_sr3_layer27, en_sr3_layer28};"
    "std::vector< std::vector<float> > en_layer = {en_layer1, en_layer2, en_layer3};"
    "return en_layer;";

  TFile *file = new TFile("file_after_weights.root", "RECREATE");
  file->cd();
  TTree *tree = new TTree("data", "Tree after weights");
  float_ geneta;
  vec1d<float_> deltaE(nreg), deltaE_corr(nreg);
  tree->Branch("geneta", &geneta);
  tree->Branch("deltaE", &deltaE);
  tree->Branch("deltaE_corr", &deltaE_corr);

  auto apply_weights = [&](float_ gen, float_ geta, 
			   vec1d<float_> en, vec1d<float_> noi,
			   vec2d<float_> en_layer, vec1d<int_> showerid,
			   vec1d<float_> f1, vec1d<float_> f2) {
    SoftwareCorrection showercorr("fileweights.root");
    vector3d<float_> weights = showercorr.weights;
    vec1d<int_> boundaries(nreg);
    vec1d<float_> lshift(nreg);
    std::string corr_mode;
    if(samples=="inner") {
      boundaries = {5, 5, 5};
      corr_mode = "left";
      lshift = {.65, .59, .48};
    }
    if(samples=="outer") {
      boundaries = {23, 23, 23};
      corr_mode = "right";
      lshift = {1., 1., 1.};
    }
    vec1d<float_> lowfact = showercorr.low_stats_factor(boundaries, corr_mode);
    for(int ireg=0; ireg<nreg; ++ireg) {
      bool_ weight_limit;
      float_ encorr = 0.;
      for (int il=1; il<=nlayers; ++il) {
	if(samples=="inner")
	  weight_limit = il > boundaries[ireg];
	else if(samples=="outer")
	  weight_limit = il < boundaries[ireg];
	float_ v = f1[ireg]*en_layer[ireg][il-1] - f2[ireg];
	if(showerid[ireg]==0) {
	  encorr += v;
	}
	else {
	  int_ w = showerid[ireg]-1;
	  if(weights(ireg,w,il-1) != 0 && weight_limit) {
	    int_ r = static_cast<int_>( std::round((il-1)*lshift[w]) );
	    encorr += v/weights(ireg,w,r);
	  }
	}
      }
      deltaE[ireg] = en[ireg]/gen - 1.;
      if(showerid[ireg]!=0) {
	encorr *= ( 1 / (1-lowfact[ireg]) );
	if(samples=="inner") encorr *= 1/0.09;
	else if(samples=="outer") encorr *= 1/0.08;
      }
      deltaE_corr[ireg] = encorr/gen - 1.;
    }
    geneta = geta;
    tree->Fill();
  };

  TFile* fa = new TFile(noPUFile.c_str(), "READ");
  TFile* fb = new TFile("filefriend.root", "READ");
  TFile* fc = new TFile("filefriend2.root", "READ");
  fa->cd();
  TTree* ta = static_cast<TTree*>(fa->Get("data"));
  fb->cd();
  TTree* tb = static_cast<TTree*>(fb->Get("tfriend1"));
  fc->cd();
  TTree* tc = static_cast<TTree*>(fc->Get("tfriend2"));
  fa->cd();
  ta->AddFriend(tb, "friend_b");
  ta->AddFriend(tc, "friend_c");

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame dfinal(*ta);
  dfinal.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(apply_weights, {"genen", "abs_geneta", "en", "noise", "en_layer", "friend_c.showerid", "friend_b.f1", "friend_b.f2"});

  file->cd();
  tree->Write();
  file->Close();
  delete file;
  return 0;
}

//use lowfact!
