f = open ("template_GenerationMixturefewersignaturespairedPCAWG.R", "r")
a = f.read()
print(a)
substitute='Stomach-AdenoCA'

f2 = open ("template_GenerationMixturefewersmallsignaturespairedPCAWG.R", "r")
a2 = f2.read()
substitute2='Kidney-RCC.clearcell'



dic = {"BoneOsteosarcPCAWG": "Bone-Osteosarc", "BreastAdenoCAPCAWG": "Breast-AdenoCA", "CNSGBMPCAWG": "CNS-GBM",\
"CNSMedulloPCAWG": "CNS-Medullo", "CNSPiloAstroPCAWG": "CNS-PiloAstro", "ColoRectAdenoCAPCAWG": "ColoRect-AdenoCA",\
"EsoAdenoCAPCAWG": "Eso-AdenoCA", "HeadSCCPCAWG": "Head-SCC", "KidneyChRCCPCAWG": "Kidney-ChRCC",\
"KidneyRCCPCAWG": "Kidney-RCC.clearcell", "KidneyRCCpapillaryPCAWG":  "Kidney-RCC.papillary",\
"LiverHCC": "Liver-HCC", "LungSCCPCAWG": "Lung-SCC", "LymphBNHLPCAWG":  "Lymph-BNHL", "LymphCLLPCAWG":  "Lymph-CLL",\
"OvaryAdenoCAPCAWG": "Ovary-AdenoCA", "PancAdenoCAPCAWG": "Panc-AdenoCA", "PancEndocrinePCAWG": "Panc-Endocrine",\
"ProstAdenoCAPCAWG": "Prost-AdenoCA", "SkinMelanomacutaneousPCAWG": "Skin-Melanoma.cutaneous",\
"StomachPCAWG": "Stomach-AdenoCA", "ThyAdenoCAPCAWG": "Thy-AdenoCA", "UterusAdenoCAPCAWG": "Uterus-AdenoCA" }

for i in dic.keys():
    o = open ("GenerationMixturefewersignaturespaired"+i+".R", "w")
    o.write(a.replace(substitute, dic[i]))
    o.close()

    o2 = open ("../../../../data/assessing_models_simulation/GenerationMixturefewersignaturespaired"+i, "w")
    o2.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
    o2.close()

    ## signatures of lowest abundance
    o3 = open ("GenerationMixturefewersmallsignaturespaired"+i+".R", "w")
    o3.write(a2.replace(substitute2, dic[i]))
    o3.close()

print(dic.keys())
