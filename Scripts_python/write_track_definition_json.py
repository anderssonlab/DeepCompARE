import json

track_info={0: "cage_hepg2_regression",
            1: "cage_k562_regression",
            2: "dhs_hepg2_regression",
            3: "dhs_k562_regression",
            4: "starr_hepg2_regression",
            5: "starr_k562_regression",
            6: "sure_hepg2_regression",
            7: "sure_k562_regression",
            8: "cage_hepg2_classification",
            9: "cage_k562_classification",
            10: "dhs_hepg2_classification",
            11: "dhs_k562_classification",
            12: "starr_hepg2_classification",
            13: "starr_k562_classification",
            14: "sure_hepg2_classification",
            15: "sure_k562_classification",
        }

# dump track_info to json file
with open('/isdata/alab/people/pcr980/DeepCompare/Scripts_python/track_definition.json', 'w') as outfile:
    json.dump(track_info, outfile)