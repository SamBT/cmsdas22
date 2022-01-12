#ifndef OUTPUT_TREE_H
#define OUTPUT_TREE_H

#include "TTree.h"

// to add a variable
// 1) add it to the list of declared member variables
// 2) add the corresponding SetBranchAddress in init()
// 3) add a clear default value in clear()

class output_tree {
    public:
        output_tree();
        ~output_tree();

        void clear_vars();

        int fill()  {return tree_->Fill();}
        int write() {return tree_->Write();}
        TTree* get_tree() {return tree_.get();}

        // here list all the variables
        unsigned int run_;
        unsigned int luminosityBlock_;
        long long    event_;
        float        xs_;

        float        btag_SF_;
        float        btag_SF_bup_;
        float        btag_SF_bdown_;
        float        btag_SF_cup_;
        float        btag_SF_cdown_;
        float        btag_SF_lightup_;
        float        btag_SF_lightdown_;
        float        norm_weight_;
        float        trigger_SF_;

        int  n_btag_;

        // Truth matched data
        bool allJets_truMatch_;

        float truMatch_H1_pt_;
        float truMatch_H1_pt_breg_;
        float truMatch_H1_eta_;
        float truMatch_H1_phi_;
        float truMatch_H1_m_;
        float truMatch_H1_m_breg_;

        float truMatch_H1_pt_gen_;
        float truMatch_H1_eta_gen_;
        float truMatch_H1_phi_gen_;
        float truMatch_H1_m_gen_;

        float truMatch_H2_pt_;
        float truMatch_H2_pt_breg_;
        float truMatch_H2_eta_;
        float truMatch_H2_phi_;
        float truMatch_H2_m_;
        float truMatch_H2_m_breg_;

        float truMatch_H2_pt_gen_;
        float truMatch_H2_eta_gen_;
        float truMatch_H2_phi_gen_;
        float truMatch_H2_m_gen_;

        float truMatch_HH_pt_;
        float truMatch_HH_pt_breg_;
        float truMatch_HH_eta_;
        float truMatch_HH_phi_;
        float truMatch_HH_m_;
        float truMatch_HH_m_breg_;

        float truMatch_HH_pt_gen_;
        float truMatch_HH_eta_gen_;
        float truMatch_HH_phi_gen_;
        float truMatch_HH_m_gen_;

        // composite candidates
        float H1_pt_;
        float H1_eta_;
        float H1_phi_;
        float H1_m_;

        float H2_pt_;
        float H2_eta_;
        float H2_phi_;
        float H2_m_;

        float HH_pt_;
        float HH_eta_;
        float HH_phi_;
        float HH_m_;

        // other variables
        float H1H2_deltaEta_;
        float H1H2_deltaPhi_;
        float H1_costhetaCM_;

        // jets
        float H1_b1_pt_;
        float H1_b1_eta_;
        float H1_b1_phi_;
        float H1_b1_m_;

        float H1_b2_pt_;
        float H1_b2_eta_;
        float H1_b2_phi_;
        float H1_b2_m_;

        float H2_b1_pt_;
        float H2_b1_eta_;
        float H2_b1_phi_;
        float H2_b1_m_;

        float H2_b2_pt_;
        float H2_b2_eta_;
        float H2_b2_phi_;
        float H2_b2_m_;

        float rndm_1_;
        float rndm_2_;
        float rndm_3_;

        // Gen info for signal
        float gen_H1_m_;
        float gen_H1_pt_;
        float gen_H1_eta_;
        float gen_H1_phi_;

        float gen_H2_m_;
        float gen_H2_pt_;
        float gen_H2_eta_;
        float gen_H2_phi_;

        float gen_mHH_;
        float gen_costh_H1_cm_;
        float gen_costh_H2_cm_;

        float gen_H1_b1_m_;
        float gen_H1_b1_pt_;
        float gen_H1_b1_eta_;
        float gen_H1_b1_phi_;

        float gen_H1_b2_m_;
        float gen_H1_b2_pt_;
        float gen_H1_b2_eta_;
        float gen_H1_b2_phi_;

        float gen_H2_b1_m_;
        float gen_H2_b1_pt_;
        float gen_H2_b1_eta_;
        float gen_H2_b1_phi_;

        float gen_H2_b2_m_;
        float gen_H2_b2_pt_;
        float gen_H2_b2_eta_;
        float gen_H2_b2_phi_;


    private:
        std::unique_ptr<TTree> tree_;
        void init();
};

output_tree::output_tree()
{
    init();
}

output_tree::~output_tree()
{
    tree_.release(); // to avoid crashes at the end of execution
}

void output_tree::init()
{
    tree_ = std::unique_ptr<TTree>(new TTree("bbbbTree", "bbbbTree"));

    // create branches

    tree_ -> Branch ("run",              &run_);
    tree_ -> Branch ("luminosityBlock",  &luminosityBlock_);
    tree_ -> Branch ("event",            &event_);
    tree_ -> Branch ("xs",               &xs_);

    tree_ ->Branch ("btag_SF",           &btag_SF_);
    tree_ ->Branch ("btag_SF_bup",       &btag_SF_bup_);
    tree_ ->Branch ("btag_SF_bdown",     &btag_SF_bdown_);
    tree_ ->Branch ("btag_SF_cup",       &btag_SF_cup_);
    tree_ ->Branch ("btag_SF_cdown",     &btag_SF_cdown_);
    tree_ ->Branch ("btag_SF_lightup",   &btag_SF_lightup_);
    tree_ ->Branch ("btag_SF_lightdown", &btag_SF_lightdown_);
    tree_ ->Branch ("norm_weight",       &norm_weight_);
    tree_ ->Branch ("trigger_SF",        &trigger_SF_);

    tree_ ->Branch ("n_btag",  &n_btag_);

    tree_ -> Branch ("allJets_truMatch", &allJets_truMatch_);

    tree_ -> Branch ("truMatch_H1_pt", &truMatch_H1_pt_);
    tree_ -> Branch ("truMatch_H1_pt_breg", &truMatch_H1_pt_breg_);
    tree_ -> Branch ("truMatch_H1_eta", &truMatch_H1_eta_);
    tree_ -> Branch ("truMatch_H1_phi", &truMatch_H1_phi_);
    tree_ -> Branch ("truMatch_H1_m", &truMatch_H1_m_);
    tree_ -> Branch ("truMatch_H1_m_breg", &truMatch_H1_m_breg_);

    tree_ -> Branch("truMatch_H1_pt_gen", &truMatch_H1_pt_gen_);
    tree_ -> Branch("truMatch_H1_eta_gen", &truMatch_H1_eta_gen_);
    tree_ -> Branch("truMatch_H1_phi_gen", &truMatch_H1_phi_gen_);
    tree_ -> Branch("truMatch_H1_m_gen", &truMatch_H1_m_gen_);

    tree_ -> Branch ("truMatch_H2_pt", &truMatch_H2_pt_);
    tree_ -> Branch ("truMatch_H2_pt_breg", &truMatch_H2_pt_breg_);
    tree_ -> Branch ("truMatch_H2_eta", &truMatch_H2_eta_);
    tree_ -> Branch ("truMatch_H2_phi", &truMatch_H2_phi_);
    tree_ -> Branch ("truMatch_H2_m", &truMatch_H2_m_);
    tree_ -> Branch ("truMatch_H2_m_breg", &truMatch_H2_m_breg_);

    tree_ -> Branch("truMatch_H2_pt_gen", &truMatch_H2_pt_gen_);
    tree_ -> Branch("truMatch_H2_eta_gen", &truMatch_H2_eta_gen_);
    tree_ -> Branch("truMatch_H2_phi_gen", &truMatch_H2_phi_gen_);
    tree_ -> Branch("truMatch_H2_m_gen", &truMatch_H2_m_gen_);

    tree_ -> Branch ("truMatch_HH_pt", &truMatch_HH_pt_);
    tree_ -> Branch ("truMatch_HH_pt_breg", &truMatch_HH_pt_breg_);
    tree_ -> Branch ("truMatch_HH_eta", &truMatch_HH_eta_);
    tree_ -> Branch ("truMatch_HH_phi", &truMatch_HH_phi_);
    tree_ -> Branch ("truMatch_HH_m", &truMatch_HH_m_);
    tree_ -> Branch ("truMatch_HH_m_breg", &truMatch_HH_m_breg_);

    tree_ -> Branch("truMatch_HH_pt_gen", &truMatch_HH_pt_gen_);
    tree_ -> Branch("truMatch_HH_eta_gen", &truMatch_HH_eta_gen_);
    tree_ -> Branch("truMatch_HH_phi_gen", &truMatch_HH_phi_gen_);
    tree_ -> Branch("truMatch_HH_m_gen", &truMatch_HH_m_gen_);

    tree_ -> Branch ("H1_pt",  &H1_pt_);
    tree_ -> Branch ("H1_eta", &H1_eta_);
    tree_ -> Branch ("H1_phi", &H1_phi_);
    tree_ -> Branch ("H1_m",   &H1_m_);

    tree_ -> Branch ("H2_pt",  &H2_pt_);
    tree_ -> Branch ("H2_eta", &H2_eta_);
    tree_ -> Branch ("H2_phi", &H2_phi_);
    tree_ -> Branch ("H2_m",   &H2_m_);

    tree_ -> Branch ("HH_pt",  &HH_pt_);
    tree_ -> Branch ("HH_eta", &HH_eta_);
    tree_ -> Branch ("HH_phi", &HH_phi_);
    tree_ -> Branch ("HH_m",   &HH_m_);

    tree_ -> Branch ("H1H2_deltaEta", &H1H2_deltaEta_);
    tree_ -> Branch ("H1H2_deltaPhi", &H1H2_deltaPhi_);
    tree_ -> Branch ("H1_costhetaCM", &H1_costhetaCM_);

    tree_ -> Branch ("H1_b1_pt",  &H1_b1_pt_);
    tree_ -> Branch ("H1_b1_eta", &H1_b1_eta_);
    tree_ -> Branch ("H1_b1_phi", &H1_b1_phi_);
    tree_ -> Branch ("H1_b1_m",   &H1_b1_m_);

    tree_ -> Branch ("H1_b2_pt",  &H1_b2_pt_);
    tree_ -> Branch ("H1_b2_eta", &H1_b2_eta_);
    tree_ -> Branch ("H1_b2_phi", &H1_b2_phi_);
    tree_ -> Branch ("H1_b2_m",   &H1_b2_m_);

    tree_ -> Branch ("H2_b1_pt",  &H2_b1_pt_);
    tree_ -> Branch ("H2_b1_eta", &H2_b1_eta_);
    tree_ -> Branch ("H2_b1_phi", &H2_b1_phi_);
    tree_ -> Branch ("H2_b1_m",   &H2_b1_m_);

    tree_ -> Branch ("H2_b2_pt",  &H2_b2_pt_);
    tree_ -> Branch ("H2_b2_eta", &H2_b2_eta_);
    tree_ -> Branch ("H2_b2_phi", &H2_b2_phi_);
    tree_ -> Branch ("H2_b2_m",   &H2_b2_m_);

    tree_ -> Branch ("rndm_1", &rndm_1_);
    tree_ -> Branch ("rndm_2", &rndm_2_);
    tree_ -> Branch ("rndm_3", &rndm_3_);

    tree_ -> Branch ("gen_H1_m", &gen_H1_m_);
    tree_ -> Branch ("gen_H1_pt", &gen_H1_pt_);
    tree_ -> Branch ("gen_H1_eta", &gen_H1_eta_);
    tree_ -> Branch ("gen_H1_phi", &gen_H1_phi_);

    tree_ -> Branch ("gen_H2_m", &gen_H2_m_);
    tree_ -> Branch ("gen_H2_pt", &gen_H2_pt_);
    tree_ -> Branch ("gen_H2_eta", &gen_H2_eta_);
    tree_ -> Branch ("gen_H2_phi", &gen_H2_phi_);

    tree_ -> Branch ("gen_mHH", &gen_mHH_);
    tree_ -> Branch ("gen_costh_H1_cm", &gen_costh_H1_cm_);
    tree_ -> Branch ("gen_costh_H2_cm", &gen_costh_H2_cm_);

    tree_ -> Branch ("gen_H1_b1_m", &gen_H1_b1_m_);
    tree_ -> Branch ("gen_H1_b1_pt", &gen_H1_b1_pt_);
    tree_ -> Branch ("gen_H1_b1_eta", &gen_H1_b1_eta_);
    tree_ -> Branch ("gen_H1_b1_phi", &gen_H1_b1_phi_);

    tree_ -> Branch ("gen_H1_b2_m", &gen_H1_b2_m_);
    tree_ -> Branch ("gen_H1_b2_pt", &gen_H1_b2_pt_);
    tree_ -> Branch ("gen_H1_b2_eta_", &gen_H1_b2_eta_);
    tree_ -> Branch ("gen_H1_b2_phi", &gen_H1_b2_phi_);

    tree_ -> Branch ("gen_H2_b1_m", &gen_H2_b1_m_);
    tree_ -> Branch ("gen_H2_b1_pt", &gen_H2_b1_pt_);
    tree_ -> Branch ("gen_H2_b1_eta", &gen_H2_b1_eta_);
    tree_ -> Branch ("gen_H2_b1_phi", &gen_H2_b1_phi_);

    tree_ -> Branch ("gen_H2_b2_m", &gen_H2_b2_m_);
    tree_ -> Branch ("gen_H2_b2_pt", &gen_H2_b2_pt_);
    tree_ -> Branch ("gen_H2_b2_eta", &gen_H2_b2_eta_);
    tree_ -> Branch ("gen_H2_b2_phi", &gen_H2_b2_phi_);

}

void output_tree::clear_vars()
{
    run_               = 0;
    luminosityBlock_   = 0;
    event_             = -999;
    xs_                = -999;

    btag_SF_           = -999;
    btag_SF_bup_       = -999;
    btag_SF_bdown_     = -999;
    btag_SF_cup_       = -999;
    btag_SF_cdown_     = -999;
    btag_SF_lightup_   = -999;
    btag_SF_lightdown_ = -999;
    norm_weight_       = -999;
    trigger_SF_        = -999;

    n_btag_ = -999;

    allJets_truMatch_ = false;

    truMatch_H1_pt_ = -999;
    truMatch_H1_pt_breg_ = -999;
    truMatch_H1_eta_ = -999;
    truMatch_H1_phi_ = -999;
    truMatch_H1_m_ = -999;
    truMatch_H1_m_breg_ = -999;
    
    truMatch_H1_pt_gen_ = -999;
    truMatch_H1_eta_gen_ = -999;
    truMatch_H1_phi_gen_ = -999;
    truMatch_H1_m_gen_ = -999;

    truMatch_H2_pt_ = -999;
    truMatch_H2_pt_breg_ = -999;
    truMatch_H2_eta_ = -999;
    truMatch_H2_phi_ = -999;
    truMatch_H2_m_ = -999;
    truMatch_H2_m_breg_ = -999;

    truMatch_H2_pt_gen_ = -999;
    truMatch_H2_eta_gen_ = -999;
    truMatch_H2_phi_gen_ = -999;
    truMatch_H2_m_gen_ = -999;

    truMatch_HH_pt_ = -999;
    truMatch_HH_pt_breg_ = -999;
    truMatch_HH_eta_ = -999;
    truMatch_HH_phi_ = -999;
    truMatch_HH_m_ = -999;
    truMatch_HH_m_breg_ = -999;

    truMatch_HH_pt_gen_ = -999;
    truMatch_HH_eta_gen_ = -999;
    truMatch_HH_phi_gen_ = -999;
    truMatch_HH_m_gen_ = -999;

    H1_pt_  = -999;
    H1_eta_ = -999;
    H1_phi_ = -999;
    H1_m_   = -999;

    H2_pt_  = -999;
    H2_eta_ = -999;
    H2_phi_ = -999;
    H2_m_   = -999;

    HH_pt_  = -999;
    HH_eta_ = -999;
    HH_phi_ = -999;
    HH_m_   = -999;

    H1H2_deltaEta_ = -999;
    H1H2_deltaPhi_ = -999;
    H1_costhetaCM_ = -999;

    H1_b1_pt_  = -999;
    H1_b1_eta_ = -999;
    H1_b1_phi_ = -999;
    H1_b1_m_   = -999;

    H1_b2_pt_  = -999;
    H1_b2_eta_ = -999;
    H1_b2_phi_ = -999;
    H1_b2_m_   = -999;

    H2_b1_pt_  = -999;
    H2_b1_eta_ = -999;
    H2_b1_phi_ = -999;
    H2_b1_m_   = -999;

    H2_b2_pt_  = -999;
    H2_b2_eta_ = -999;
    H2_b2_phi_ = -999;
    H2_b2_m_   = -999;

    rndm_1_ = -999;
    rndm_2_ = -999;
    rndm_3_ = -999;

    gen_H1_m_ = -999;
    gen_H1_pt_ = -999;
    gen_H1_eta_ = -999;
    gen_H1_phi_ = -999;

    gen_H2_m_ = -999;
    gen_H2_pt_ = -999;
    gen_H2_eta_ = -999;
    gen_H2_phi_ = -999;

    gen_mHH_ = -999;
    gen_costh_H1_cm_ = -999;
    gen_costh_H2_cm_ = -999;

    gen_H1_b1_m_ = -999;
    gen_H1_b1_pt_ = -999;
    gen_H1_b1_eta_ = -999;
    gen_H1_b1_phi_ = -999;

    gen_H1_b2_m_ = -999;
    gen_H1_b2_pt_ = -999;
    gen_H1_b2_eta_ = -999;
    gen_H1_b2_phi_ = -999;

    gen_H2_b1_m_ = -999;
    gen_H2_b1_pt_ = -999;
    gen_H2_b1_eta_ = -999;
    gen_H2_b1_phi_ = -999;

    gen_H2_b2_m_ = -999;
    gen_H2_b2_pt_ = -999;
    gen_H2_b2_eta_ = -999;
    gen_H2_b2_phi_ = -999;
}

#endif