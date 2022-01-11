// puts together some properties of a jet as stored in the input file
// in this way operations like sorting etc. are easier

#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "TLorentzVector.h"

struct jet_t {

    float pt;
    float eta;
    float phi;
    float m;

    float btagscore;
    float bRegCorr;
    float bRegRes;
    int   jetId;
    int   puId;

    float genjet_pt;
    float genjet_eta;
    float genjet_phi;
    float genjet_m;
    int   genjet_hadronFlavour;

    // extra things that can be useful. NB: not initiated in the init_ function
    float pt_breg;
    float m_breg;
    TLorentzVector p4;
    TLorentzVector p4_breg;
};

// simple preprocessing macro to initialise the content from the tree
// just avoids copying stuff by hand
#define init_jet_t(OBJ, IDX, TREEREADER, SAMPLETYPE) \
    OBJ . pt                   = **( TREEREADER .  jet ## IDX ## _pt                    ) ; \
    OBJ . eta                  = **( TREEREADER .  jet ## IDX ## _eta                   ) ; \
    OBJ . phi                  = **( TREEREADER .  jet ## IDX ## _phi                   ) ; \
    OBJ . m                    = **( TREEREADER .  jet ## IDX ## _m                     ) ; \
    OBJ . btagscore            = **( TREEREADER .  jet ## IDX ## _btagscore             ) ; \
    OBJ . bRegCorr             = **( TREEREADER .  jet ## IDX ## _bRegCorr              ) ; \
    OBJ . bRegRes              = **( TREEREADER .  jet ## IDX ## _bRegRes               ) ; \
    OBJ . jetId                = **( TREEREADER .  jet ## IDX ## _jetId                 ) ; \
    OBJ . puId                 = **( TREEREADER .  jet ## IDX ## _puId                  ) ; \
    OBJ . genjet_pt            = ( SAMPLETYPE == "data" ? -999 : **( TREEREADER .  genjet ## IDX ## _pt                 ) ); \
    OBJ . genjet_eta           = ( SAMPLETYPE == "data" ? -999 : **( TREEREADER .  genjet ## IDX ## _eta                ) ); \
    OBJ . genjet_phi           = ( SAMPLETYPE == "data" ? -999 : **( TREEREADER .  genjet ## IDX ## _phi                ) ); \
    OBJ . genjet_m             = ( SAMPLETYPE == "data" ? -999 : **( TREEREADER .  genjet ## IDX ## _m                  ) ); \
    OBJ . genjet_hadronFlavour = ( SAMPLETYPE == "data" ? -999 : **( TREEREADER .  genjet ## IDX ## _hadronFlavour      ) ); \

#define init_genH_t(OBJ, IDX, TREEREADER) \
    OBJ . pt                   = **( TREEREADER .  gen_H ## IDX ## _pt                    ) ; \
    OBJ . eta                  = **( TREEREADER .  gen_H ## IDX ## _eta                   ) ; \
    OBJ . phi                  = **( TREEREADER .  gen_H ## IDX ## _phi                   ) ; \
    OBJ . m                    = **( TREEREADER .  gen_H ## IDX ## _m                     ) ; \

#define init_genHb_t(OBJ, IDXH, IDXb, TREEREADER) \
    OBJ . pt                   = **( TREEREADER .  gen_H ## IDXH ## _b ## IDXb ## _pt                    ) ; \
    OBJ . eta                  = **( TREEREADER .  gen_H ## IDXH ## _b ## IDXb ## _eta                   ) ; \
    OBJ . phi                  = **( TREEREADER .  gen_H ## IDXH ## _b ## IDXb ## _phi                   ) ; \
    OBJ . m                    = **( TREEREADER .  gen_H ## IDXH ## _b ## IDXb ## _m                     ) ; \

#endif

