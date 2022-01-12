// code used to build the objects starting from the basic inputs of the analysis

// c++ -lm -o build_objects build_objects.cpp -I include/ `root-config --glibs --cflags`

#include<iostream>
#include<string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"

#include "analysis_utils.h" // jet_t, init_jet_t
#include "input_tree.h"
#include "output_tree.h"
#include "bbbb_functions.h"

using namespace std;

int main(int argc, char** argv)
{
    // handling the I/O

    if (argc < 3)
    {
        cout << "*** usage : ./build_objects inputFile outputFile isData=0 isSignal=0" << endl;
        return 1;
    }

    string inputFile  = argv[1];
    string outputFile = argv[2];

    cout << "[INFO] Input  file: " << inputFile << endl;
    TFile* fIn = TFile::Open(inputFile.c_str());
    TTree* tIn = (TTree*) fIn->Get("bbbbTree");

    TriggerEfficiencyCalculator theTriggerEfficiencyCalculator("../trigger/TriggerEfficiencies.root");

    cout << "[INFO] Output file: " << outputFile << endl;
    TFile* fOut = new TFile(outputFile.c_str(), "recreate");
    TTree* tOut = new TTree ("bbbbTree", "bbbbTree");

    bool isData = false;
    bool isSig  = false;
    string mode = "minMassDist";
    if (argc >= 4)
        isData = (std::stoi(argv[3]) == 0 ? false : true);
    if (argc >= 5)
        isSig = (std::stoi(argv[4]) == 0 ? false : true);
    if (argc >= 6)
        mode = argv[5];
    cout << "[INFO] Is data?   : " << std::boolalpha << isData << std::noboolalpha << endl;
    cout << "[INFO] Is signal? : " << std::boolalpha << isSig  << std::noboolalpha << endl;

    std::string sample_type = "bkg";
    if (isData) sample_type = "data";
    if (isSig)  sample_type = "sig";
    cout << "[INFO] The tree will be read for type : " << sample_type << endl;

    const float btag_WP_medium = 0.3093; // for DeepJet

    // declaring the output tree - use the output_tree class
    // output variables are defined inside there
    output_tree otree;

    // prepare the input tree reader and run the event loop
    input_tree itree(tIn, sample_type);

    const auto nEv = tIn->GetEntries();
    for (uint iEv = 0; iEv < nEv; ++iEv)
    {
        tIn->GetEntry(iEv);

        if (iEv % 10000 == 0)
            cout << "... doing event " << iEv << " / " << nEv << endl;

        otree.clear_vars();

        // prepare the objects of the 4 jets
        jet_t jet1, jet2, jet3, jet4;
        TLorentzVector H1, H2, H1b1, H1b2, H2b1, H2b2;
        if (isSig) {
            H1.SetPtEtaPhiM(**(itree.gen_H1_pt),**(itree.gen_H1_eta),**(itree.gen_H1_phi),**(itree.gen_H1_m));
            H2.SetPtEtaPhiM(**(itree.gen_H2_pt),**(itree.gen_H2_eta),**(itree.gen_H2_phi),**(itree.gen_H2_m));
            H1b1.SetPtEtaPhiM(**(itree.gen_H1_b1_pt),**(itree.gen_H1_b1_eta),**(itree.gen_H1_b1_phi),**(itree.gen_H1_b1_m));
            H1b2.SetPtEtaPhiM(**(itree.gen_H1_b2_pt),**(itree.gen_H1_b2_eta),**(itree.gen_H1_b2_phi),**(itree.gen_H1_b2_m));
            H2b1.SetPtEtaPhiM(**(itree.gen_H2_b1_pt),**(itree.gen_H2_b1_eta),**(itree.gen_H2_b1_phi),**(itree.gen_H2_b1_m));
            H2b2.SetPtEtaPhiM(**(itree.gen_H2_b2_pt),**(itree.gen_H2_b2_eta),**(itree.gen_H2_b2_phi),**(itree.gen_H2_b2_m));
        }

        // the macros below copy the properties of the jet with the specified index (1,2,3,4) to the jet object (jet1,jet2,jet3,jet4)
        init_jet_t(jet1, 1, itree, sample_type) ;
        init_jet_t(jet2, 2, itree, sample_type) ;
        init_jet_t(jet3, 3, itree, sample_type) ;
        init_jet_t(jet4, 4, itree, sample_type) ;


        // compute regressed quantities and other member variables
        jet1.pt_breg = jet1.bRegCorr * jet1.pt;
        jet1.m_breg  = jet1.bRegCorr * jet1.m;
        jet1.p4      .SetPtEtaPhiM(jet1.pt,      jet1.eta, jet1.phi, jet1.m);
        jet1.p4_breg .SetPtEtaPhiM(jet1.pt_breg, jet1.eta, jet1.phi, jet1.m_breg);

        jet2.pt_breg = jet2.bRegCorr * jet2.pt;
        jet2.m_breg  = jet2.bRegCorr * jet2.m;
        jet2.p4      .SetPtEtaPhiM(jet2.pt,      jet2.eta, jet2.phi, jet2.m);
        jet2.p4_breg .SetPtEtaPhiM(jet2.pt_breg, jet2.eta, jet2.phi, jet2.m_breg);

        jet3.pt_breg = jet3.bRegCorr * jet3.pt;
        jet3.m_breg  = jet3.bRegCorr * jet3.m;
        jet3.p4      .SetPtEtaPhiM(jet3.pt,      jet3.eta, jet3.phi, jet3.m);
        jet3.p4_breg .SetPtEtaPhiM(jet3.pt_breg, jet3.eta, jet3.phi, jet3.m_breg);

        jet4.pt_breg = jet4.bRegCorr * jet4.pt;
        jet4.m_breg  = jet4.bRegCorr * jet4.m;
        jet4.p4      .SetPtEtaPhiM(jet4.pt,      jet4.eta, jet4.phi, jet4.m);
        jet4.p4_breg .SetPtEtaPhiM(jet4.pt_breg, jet4.eta, jet4.phi, jet4.m_breg);

        // ========================================
        // ========================================
        // ---- here goes your analysis code to build the high level objects
        // cout << "JET 1 : " << jet1.pt << " " << jet1.eta << " " << jet1.genjet_pt << " " << endl;
        // cout << "JET 2 : " << jet2.pt << " " << jet2.eta << " " << jet2.genjet_pt << " " << endl;
        // cout << "JET 3 : " << jet3.pt << " " << jet3.eta << " " << jet3.genjet_pt << " " << endl;
        // cout << "JET 4 : " << jet4.pt << " " << jet4.eta << " " << jet4.genjet_pt << " " << endl;
        // cout << " -------------------------------------- " << endl;

        // For Monte Carlo, match jets to gen-level b-quarks
        if (isSig) {
            float Rt = 0.1;
            vector<TLorentzVector> bq = {H1b1,H1b2,H2b1,H2b2};
            vector<jet_t> js = {jet1,jet2,jet3,jet4};
            std::map<int,jet_t> bquark_jMatch;
            vector<int> matched_bqs;
            jet_t jcurr;
            TLorentzVector pcurr;
            int n_match = 0;
            for (int ij = 0; ij < 4; ij++) {
                jcurr = js[ij];
                if (jcurr.genjet_hadronFlavour != 5) continue;
                float dRmin = 999;
                int imatch = -1;
                for (int ip = 0; ip < 4; ip++) {
                    pcurr = bq[ip];
                    float dr = jcurr.p4.DeltaR(pcurr);
                    if (dr < dRmin) {
                        dRmin = dr;
                        imatch = ip;
                    }
                }
                if (imatch != -1 && dRmin < Rt && std::find(matched_bqs.begin(),matched_bqs.end(),imatch) == matched_bqs.end()) {
                    bquark_jMatch[imatch] = js[ij];
                    matched_bqs.push_back(imatch);
                    n_match++;
                }
            }
            if (n_match == 4) {
                otree.allJets_truMatch_ = true;
                TLorentzVector true_H1 = bquark_jMatch[0].p4 + bquark_jMatch[1].p4;
                TLorentzVector true_H2 = bquark_jMatch[2].p4 + bquark_jMatch[3].p4;
                TLorentzVector true_HH = true_H1 + true_H2;

                TLorentzVector true_H1_breg = bquark_jMatch[0].p4_breg + bquark_jMatch[1].p4_breg;
                TLorentzVector true_H2_breg = bquark_jMatch[2].p4_breg + bquark_jMatch[3].p4_breg;
                TLorentzVector true_HH_breg = true_H1 + true_H2;

                if (true_H1.Pt() < true_H2.Pt()) {
                    std::swap(true_H1,true_H2);
                }
                if (true_H1_breg.Pt() < true_H2_breg.Pt()) {
                    std::swap(true_H1_breg,true_H2_breg);
                }

                otree.truMatch_H1_pt_ = true_H1.Pt();
                otree.truMatch_H1_pt_breg_ = true_H1_breg.Pt();
                otree.truMatch_H1_eta_ = true_H1.Eta();
                otree.truMatch_H1_phi_ = true_H1.Phi();
                otree.truMatch_H1_m_ = true_H1.M();
                otree.truMatch_H1_m_breg_ = true_H1_breg.M();

                otree.truMatch_H2_pt_ = true_H2.Pt();
                otree.truMatch_H2_pt_breg_ = true_H2_breg.Pt();
                otree.truMatch_H2_eta_ = true_H2.Eta();
                otree.truMatch_H2_phi_ = true_H2.Phi();
                otree.truMatch_H2_m_ = true_H2.M();
                otree.truMatch_H2_m_breg_ = true_H2_breg.M();
                
                otree.truMatch_HH_pt_ = true_HH.Pt();
                otree.truMatch_HH_pt_breg_ = true_HH_breg.Pt();
                otree.truMatch_HH_eta_ = true_HH.Eta();
                otree.truMatch_HH_phi_ = true_HH.Phi();
                otree.truMatch_HH_m_ = true_HH.M();
                otree.truMatch_HH_m_breg_ = true_HH_breg.M();
            }
        }
        
        // pair the jets
        std::vector<jet_t> jets {jet1, jet2, jet3, jet4};
        std::vector<jet_t> result = bbbb_pairing(&jets,mode);

        TLorentzVector v_H1, v_H2, v_HH;
        v_H1 = result.at(0).p4_breg + result.at(1).p4_breg;
        v_H2 = result.at(2).p4_breg + result.at(3).p4_breg;

        // order them by highest pT: pt(H1) > pt(H2)

        jet_t H1_b1 = result.at(0);
        jet_t H1_b2 = result.at(1);
        jet_t H2_b1 = result.at(2);
        jet_t H2_b2 = result.at(3);

        if (v_H1.Pt() < v_H2.Pt()){
            std::swap(v_H1, v_H2);
            std::swap(H1_b1, H2_b1);
            std::swap(H1_b2, H2_b2);
        }

        v_HH = v_H1 + v_H2;
        // ========================================
        // ========================================

        // copy to the output
        otree.H1_pt_  = v_H1.Pt();
        otree.H1_eta_ = v_H1.Eta();
        otree.H1_phi_ = v_H1.Phi();
        otree.H1_m_   = v_H1.M();

        otree.H2_pt_  = v_H2.Pt();
        otree.H2_eta_ = v_H2.Eta();
        otree.H2_phi_ = v_H2.Phi();
        otree.H2_m_   = v_H2.M();

        otree.HH_pt_  = v_HH.Pt();
        otree.HH_eta_ = v_HH.Eta();
        otree.HH_phi_ = v_HH.Phi();
        otree.HH_m_   = v_HH.M();

        otree.H1_b1_pt_  = H1_b1.p4_breg.Pt();
        otree.H1_b1_eta_ = H1_b1.p4_breg.Eta();
        otree.H1_b1_phi_ = H1_b1.p4_breg.Phi();
        otree.H1_b1_m_   = H1_b1.p4_breg.M();

        otree.H1_b2_pt_  = H1_b2.p4_breg.Pt();
        otree.H1_b2_eta_ = H1_b2.p4_breg.Eta();
        otree.H1_b2_phi_ = H1_b2.p4_breg.Phi();
        otree.H1_b2_m_   = H1_b2.p4_breg.M();

        otree.H2_b1_pt_  = H2_b1.p4_breg.Pt();
        otree.H2_b1_eta_ = H2_b1.p4_breg.Eta();
        otree.H2_b1_phi_ = H2_b1.p4_breg.Phi();
        otree.H2_b1_m_   = H2_b1.p4_breg.M();

        otree.H2_b2_pt_  = H2_b2.p4_breg.Pt();
        otree.H2_b2_eta_ = H2_b2.p4_breg.Eta();
        otree.H2_b2_phi_ = H2_b2.p4_breg.Phi();
        otree.H2_b2_m_   = H2_b2.p4_breg.M();

        otree.H1H2_deltaEta_   = std::abs( v_H1.Eta() - v_H2.Eta() );
        otree.H1H2_deltaPhi_   = v_H1.DeltaPhi(v_H2);

        // boost H1 to the bbbb CM
        TLorentzVector vH1_cm = v_H1;
        vH1_cm.Boost(-v_HH.BoostVector());
        otree.H1_costhetaCM_  = vH1_cm.CosTheta();

        // and fw from the input
        otree.run_                = **(itree.run);
        otree.luminosityBlock_    = **(itree.luminosityBlock);
        otree.event_              = **(itree.event);
        otree.xs_                 = **(itree.xs);

        otree.btag_SF_            = **(itree.btag_SF);
        otree.btag_SF_bup_        = **(itree.btag_SF_bup);
        otree.btag_SF_bdown_      = **(itree.btag_SF_bdown);
        otree.btag_SF_cup_        = **(itree.btag_SF_cup);
        otree.btag_SF_cdown_      = **(itree.btag_SF_cdown);
        otree.btag_SF_lightup_    = **(itree.btag_SF_lightup);
        otree.btag_SF_lightdown_  = **(itree.btag_SF_lightdown);
        otree.norm_weight_        = **(itree.norm_weight);

        otree.n_btag_ = 0;
        if (result.at(0).btagscore > btag_WP_medium) ++otree.n_btag_;
        if (result.at(1).btagscore > btag_WP_medium) ++otree.n_btag_;
        if (result.at(2).btagscore > btag_WP_medium) ++otree.n_btag_;
        if (result.at(3).btagscore > btag_WP_medium) ++otree.n_btag_;

        otree.rndm_1_ = **(itree.rndm_1);
        otree.rndm_2_ = **(itree.rndm_2);
        otree.rndm_3_ = **(itree.rndm_3);

        if (isSig) {
            otree.gen_H1_m_ = **(itree.gen_H1_m);
            otree.gen_H1_pt_ = **(itree.gen_H1_pt);
            otree.gen_H1_eta_ = **(itree.gen_H1_eta);
            otree.gen_H1_phi_ = **(itree.gen_H1_phi);

            otree.gen_H2_m_ = **(itree.gen_H2_m);
            otree.gen_H2_pt_ = **(itree.gen_H2_pt);
            otree.gen_H2_eta_ = **(itree.gen_H2_eta);
            otree.gen_H2_phi_ = **(itree.gen_H2_phi);

            otree.gen_mHH_ = **(itree.gen_mHH);
            otree.gen_costh_H1_cm_ = **(itree.gen_costh_H1_cm);
            otree.gen_costh_H2_cm_ = **(itree.gen_costh_H2_cm);

            otree.gen_H1_b1_m_ = **(itree.gen_H1_b1_m);
            otree.gen_H1_b1_pt_ = **(itree.gen_H1_b1_pt);
            otree.gen_H1_b1_eta_ = **(itree.gen_H1_b1_eta);
            otree.gen_H1_b1_phi_ = **(itree.gen_H1_b1_phi);

            otree.gen_H1_b2_m_ = **(itree.gen_H1_b2_m);
            otree.gen_H1_b2_pt_ = **(itree.gen_H1_b2_pt);
            otree.gen_H1_b2_eta_ = **(itree.gen_H1_b2_eta);
            otree.gen_H1_b2_phi_ = **(itree.gen_H1_b2_phi);

            otree.gen_H2_b1_m_ = **(itree.gen_H2_b1_m);
            otree.gen_H2_b1_pt_ = **(itree.gen_H2_b1_pt);
            otree.gen_H2_b1_eta_ = **(itree.gen_H2_b1_eta);
            otree.gen_H2_b1_phi_ = **(itree.gen_H2_b1_phi);

            otree.gen_H2_b2_m_ = **(itree.gen_H2_b2_m);
            otree.gen_H2_b2_pt_ = **(itree.gen_H2_b2_pt);
            otree.gen_H2_b2_eta_ = **(itree.gen_H2_b2_eta);
            otree.gen_H2_b2_phi_ = **(itree.gen_H2_b2_phi);
        }
        
        // calculate the triggerSF following the twiki indications
        if(isData) otree.trigger_SF_ = 1.;
        else
        {
            // Pay attention to provide the correct variable to estract the trigger efficiency!!
            std::vector<float> jetPtVector {jet1.pt, jet2.pt, jet3.pt, jet4.pt};
            // Remember to order the 4 jets by pT!! (search for std::sort)
            // std::sort...; 

            // Estract the efficiency for the four filters considered in data
            float dataEfficiency_Double90Quad30_QuadCentralJet30   = theTriggerEfficiencyCalculator.getDataEfficiency_Double90Quad30_QuadCentralJet30  (-999.);
            float dataEfficiency_Double90Quad30_DoubleCentralJet90 = theTriggerEfficiencyCalculator.getDataEfficiency_Double90Quad30_DoubleCentralJet90(-999.);
            float dataEfficiency_Quad45_QuadCentralJet45           = theTriggerEfficiencyCalculator.getDataEfficiency_Quad45_QuadCentralJet45          (-999.);
            float dataEfficiency_And_QuadCentralJet45              = theTriggerEfficiencyCalculator.getDataEfficiency_And_QuadCentralJet45             (-999.);
            // Calculate data total efficiency
            float dataEfficiency_Double90Quad30 = -999.;
            float dataEfficiency = -999.;

            // Estract the efficiency for the four filters considered in mc
            float mcEfficiency_Double90Quad30_QuadCentralJet30     = theTriggerEfficiencyCalculator.getMcEfficiency_Double90Quad30_QuadCentralJet30    (-999.);
            float mcEfficiency_Double90Quad30_DoubleCentralJet90   = theTriggerEfficiencyCalculator.getMcEfficiency_Double90Quad30_DoubleCentralJet90  (-999.);
            float mcEfficiency_Quad45_QuadCentralJet45             = theTriggerEfficiencyCalculator.getMcEfficiency_Quad45_QuadCentralJet45            (-999.);
            float mcEfficiency_And_QuadCentralJet45                = theTriggerEfficiencyCalculator.getMcEfficiency_And_QuadCentralJet45               (-999.);
            // Calculate data total efficiency
            float mcEfficiency_Double90Quad30 = -999.;
            float mcEfficiency = -999.;

            // Calculate the trigger scale factor (data/mc)
            otree.trigger_SF_ = -999.;
        }
        

        otree.fill();
    }

    // save to output
    fOut->cd();
    otree.write();

    fIn->Close();
    fOut->Close();

}