"""
Colorectal Machine Learning Final Pipeline

Instructions:
  Pipeline:
    load_and_pickle| uses cr_columns (cuts down columns) & pickles main dataframe
    pickle_comp| reads in main df from load an pickle &
    pickle_surgeries| uses main df from load and pickle & pickles surgeries

@author: Robert Handzel
Last Modified: 4/1/17 
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from matplotlib import style
style.use('ggplot')
import re

def running_fxn(splits,percent):
        print('0%|'+'#'*int(percent/(100/splits))+' '*int((100-percent)/(100/splits))+'|100%')

def load_and_pickle(file):
    print('load_and_pickle function is running...')
    
    #load and pickle cr database
    # df = pd.read_excel(file,sheetname='CR_all')
    # pd.to_pickle(df, 'S:\ERAS\cr_datebase.pickle')

    #define and pickle sx_list
    sx_list = ['redcap_event_name','patient_id','prim_sx_rectalca_a___7','prim_sx_rectalca_a___8','prim_sx_rectalca_a___9','prim_sx_rectalca_a___10','prim_sx_rectalca_a___25','prim_sx_rectalca_a___11','prim_sx_rectalca_a___31','prim_sx_rectalca_a___30','prim_sx_rectalca_a___13','prim_sx_rectalca_a___14','prim_sx_rectalca_a___15','prim_sx_rectalca_a___27','prim_sx_rectalca_a___24','prim_sx_rectalca_a___16','prim_sx_rectalca_a___17','prim_sx_rectalca_a___18','prim_sx_rectalca_a___19','prim_sx_rectalca_a___28','prim_sx_rectalca_a___29','prim_sx_rectalca_a___20','prim_sx_rectalca_a___21','prim_sx_rectalca_a___22','prim_sx_rectalca_a___23','prim_sx_other_rectalca_a',
    'prim_sx_rectalpolyp_a___7','prim_sx_rectalpolyp_a___8','prim_sx_rectalpolyp_a___9','prim_sx_rectalpolyp_a___10','prim_sx_rectalpolyp_a___25','prim_sx_rectalpolyp_a___11','prim_sx_rectalpolyp_a___12','prim_sx_rectalpolyp_a___30','prim_sx_rectalpolyp_a___29','prim_sx_rectalpolyp_a___13','prim_sx_rectalpolyp_a___26','prim_sx_rectalpolyp_a___14','prim_sx_rectalpolyp_a___15','prim_sx_rectalpolyp_a___16','prim_sx_rectalpolyp_a___24','prim_sx_rectalpolyp_a___27','prim_sx_rectalpolyp_a___17','prim_sx_rectalpolyp_a___18','prim_sx_rectalpolyp_a___19','prim_sx_rectalpolyp_a___28','prim_sx_rectalpolyp_a___20','prim_sx_rectalpolyp_a___21','prim_sx_rectalpolyp_a___22','prim_sx_rectalpolyp_a___23','prim_sx_other_rectlpolyp_a',
    'prim_sx_colonca_a___7','prim_sx_colonca_a___8','prim_sx_colonca_a___9','prim_sx_colonca_a___10','prim_sx_colonca_a___11','prim_sx_colonca_a___12','prim_sx_colonca_a___32','prim_sx_colonca_a___13','prim_sx_colonca_a___14','prim_sx_colonca_a___15','prim_sx_colonca_a___16','prim_sx_colonca_a___35','prim_sx_colonca_a___36','prim_sx_colonca_a___34','prim_sx_colonca_a___29','prim_sx_colonca_a___28','prim_sx_colonca_a___17','prim_sx_colonca_a___18','prim_sx_colonca_a___19','prim_sx_colonca_a___27','prim_sx_colonca_a___20','prim_sx_colonca_a___30','prim_sx_colonca_a___21','prim_sx_colonca_a___22','prim_sx_colonca_a___31','prim_sx_colonca_a___23','prim_sx_colonca_a___24','prim_sx_colonca_a___25','prim_sx_colonca_a___26','prim_sx_other_colonca_a',
    'prim_sx_colonpolyp_a___7','prim_sx_colonpolyp_a___8','prim_sx_colonpolyp_a___9','prim_sx_colonpolyp_a___10','prim_sx_colonpolyp_a___11','prim_sx_colonpolyp_a___12','prim_sx_colonpolyp_a___32','prim_sx_colonpolyp_a___13','prim_sx_colonpolyp_a___14','prim_sx_colonpolyp_a___15','prim_sx_colonpolyp_a___16','prim_sx_colonpolyp_a___33','prim_sx_colonpolyp_a___34','prim_sx_colonpolyp_a___35','prim_sx_colonpolyp_a___29','prim_sx_colonpolyp_a___28','prim_sx_colonpolyp_a___17','prim_sx_colonpolyp_a___18','prim_sx_colonpolyp_a___19','prim_sx_colonpolyp_a___20','prim_sx_colonpolyp_a___30','prim_sx_colonpolyp_a___21','prim_sx_colonpolyp_a___27','prim_sx_colonpolyp_a___22','prim_sx_colonpolyp_a___31','prim_sx_colonpolyp_a___23','prim_sx_colonpolyp_a___24','prim_sx_colonpolyp_a___25','prim_sx_colonpolyp_a___26','prim_sx_other_colonpolyp_a',
    'prim_sx_bencolon_a___1','prim_sx_bencolon_a___2','prim_sx_bencolon_a___3','prim_sx_bencolon_a___4','prim_sx_bencolon_a___5','prim_sx_bencolon_a___6','prim_sx_bencolon_a___7','prim_sx_bencolon_a___8','prim_sx_bencolon_a___9','prim_sx_bencolon_a___10','prim_sx_bencolon_a___11','prim_sx_bencolon_a___12','prim_sx_bencolon_a___13','prim_sx_bencolon_a___14','prim_sx_bencolon_a___15','prim_sx_bencolon_a___16','prim_sx_bencolon_a___25','prim_sx_bencolon_a___17','prim_sx_bencolon_a___18','prim_sx_bencolon_a___24','prim_sx_bencolon_a___19','prim_sx_bencolon_a___26','prim_sx_bencolon_a___27','prim_sx_bencolon_a___20','prim_sx_bencolon_a___21','prim_sx_bencolon_a___22','prim_sx_bencolon_a___23','prim_sx_other_bencolon_a',
    'sx_rectopexy_a',
    'prim_sx_uc_a___1','prim_sx_uc_a___2','prim_sx_uc_a___31','prim_sx_uc_a___3','prim_sx_uc_a___4','prim_sx_uc_a___5','prim_sx_uc_a___24','prim_sx_uc_a___25','prim_sx_uc_a___6','prim_sx_uc_a___7','prim_sx_uc_a___8','prim_sx_uc_a___29','prim_sx_uc_a___9','prim_sx_uc_a___10','prim_sx_uc_a___11','prim_sx_uc_a___22','prim_sx_uc_a___23','prim_sx_uc_a___12','prim_sx_uc_a___13','prim_sx_uc_a___14','prim_sx_uc_a___15','prim_sx_uc_a___26','prim_sx_uc_a___16','prim_sx_uc_a___21','prim_sx_uc_a___27','prim_sx_uc_a___28','prim_sx_uc_a___17','prim_sx_uc_a___18','prim_sx_uc_a___19','prim_sx_uc_a___20','prim_sx_other_uc_a',
    'prim_sx_ic_a___1','prim_sx_ic_a___2','prim_sx_ic_a___30','prim_sx_ic_a___31','prim_sx_ic_a___3','prim_sx_ic_a___4','prim_sx_ic_a___24','prim_sx_ic_a___25','prim_sx_ic_a___5','prim_sx_ic_a___6','prim_sx_ic_a___7','prim_sx_ic_a___8','prim_sx_ic_a___29','prim_sx_ic_a___9','prim_sx_ic_a___10','prim_sx_ic_a___11','prim_sx_ic_a___22','prim_sx_ic_a___23','prim_sx_ic_a___12','prim_sx_ic_a___13','prim_sx_ic_a___14','prim_sx_ic_a___15','prim_sx_ic_a___26','prim_sx_ic_a___16','prim_sx_ic_a___21','prim_sx_ic_a___27','prim_sx_ic_a___28','prim_sx_ic_a___17','prim_sx_ic_a___18','prim_sx_ic_a___19','prim_sx_ic_a___20','prim_sx_other_ic_a',
    'prim_sx_cd_a___1','prim_sx_cd_a___2','prim_sx_cd_a___30','prim_sx_cd_a___31','prim_sx_cd_a___3','prim_sx_cd_a___4','prim_sx_cd_a___24','prim_sx_cd_a___25','prim_sx_cd_a___5','prim_sx_cd_a___6','prim_sx_cd_a___7','prim_sx_cd_a___8','prim_sx_cd_a___29','prim_sx_cd_a___9','prim_sx_cd_a___10','prim_sx_cd_a___11','prim_sx_cd_a___22','prim_sx_cd_a___23','prim_sx_cd_a___26','prim_sx_cd_a___12','prim_sx_cd_a___13','prim_sx_cd_a___14','prim_sx_cd_a___15','prim_sx_cd_a___16','prim_sx_cd_a___21','prim_sx_cd_a___27','prim_sx_cd_a___28','prim_sx_cd_a___17','prim_sx_cd_a___18','prim_sx_cd_a___19','prim_sx_cd_a___20','prim_sx_other_cd_a',
    'sx_multivisc_rxn_a','sx_anastomosis_a','sx_anastamosis_ibd_a','sx_temp_diversion_a','secondary_sx_a___17','secondary_sx_a___18','secondary_sx_a___19','secondary_sx_a___20','secondary_sx_a___21','secondary_sx_a___22','secondary_sx_a___23','secondary_sx_a___24','secondary_sx_a___25','secondary_sx_a___26','secondary_sx_a___27','secondary_sx_a___28','secondary_sx_a___29','secondary_sx_a___30','other_secondary_sx_a']

    pd.to_pickle(sx_list, 'S:\ERAS\crdb_sx_list.pickle')

def extract_sx_data():
    df = pd.read_pickle('S:\ERAS\cr_datebase.pickle')
    sx_list = list(pd.read_pickle('S:\ERAS\crdb_sx_list.pickle'))
    
    redcap_events = ['surgery_dx_1_arm_1']
    sx_df = df[sx_list]
    sx_df = sx_df[sx_df.redcap_event_name.isin(redcap_events)]
    sx_df.drop('redcap_event_name',axis=1,inplace=True) #drops redcap column
    pd.to_pickle(sx_df, 'S:\ERAS\crdb_sx_data.pickle')

def create_sx_values():
    print('create_sx_values is running...')
    df = pd.read_pickle('S:\ERAS\crdb_sx_data.pickle')
    sx_ddf = pd.read_pickle('S:\ERAS\sx_list_dict_comp.pickle') #dictionary dataframe
    pt_id = list(df.patient_id)

    cnt_list = []
    num_of_pts = df.patient_id.unique().shape[0]

    sx_list = list(sx_ddf.unique.unique())

    #loops through all unique surgeries and adds a column
    #will function like a onehotencoder. each surgery will have 1 column and 0 or 1
    #one patient had 7 surgeries listed
    # for sx in sx_list:
    #     df[sx] = np.NaN
        #need to replace characters to allow for lists
        # sx = sx.replace(' ','_')
        # sx = sx.replace('/','_')
        # sx = sx.replace('-','_')
        
        # print("'{}':[],".format(sx)) #used to create text to initiate all lists
        # print("{}".format(sx)) #used to create text to initiate all lists
    
    
    output_dict = {'patient_id':[],
    'abdominal_perineal_resection':[],
    'low_anterior_resection':[],
    'proctectomy':[],
    'proctosigmoidectomy':[],
    'proctosigmoidectomy_with_colonic_j_pouch_anal_anastomosis':[],
    'transanal_endoscopic_microsurgery_transanal_minimally_invasive_surgery':[],
    'total_proctocolectomy':[],
    'total_proctocolectomy_with_end_ileostomy':[],
    'transanal_excision':[],
    'colostomy_reversal':[],
    'ileostomy_reversal':[],
    'exploratory_laparotomy':[],
    'anastomotic_leak_repair':[],
    'diverting_stoma':[],
    'incisional_hernia_repair':[],
    'parastomal_hernia_repair':[],
    'lysis_of_adhesions':[],
    'omentectomy':[],
    'enterotomy':[],
    'treatment_of_intrabdominal_bleeding':[],
    'wound_dehiscence_repair':[],
    'wound_exploration_debridement':[],
    'other':[],
    'left_hemi_colectomy':[],
    'right_hemi_colectomy':[],
    'colectomy':[],
    'sigmoid_colectomy':[],
    'transverse_colectomy':[],
    'total_abdominal_colectomy':[],
    'total_abdominal_colectomy_with_ileorectal_anastamosis':[],
    'ileocolonic_resection':[],
    'altemeier':[],
    'delorme':[],
    'hartmann_procedure':[],
    'lavage_and_drain_placement':[],
    'proctocolectomy':[],
    'rectopexy':[],
    'resection_with_rectopexy':[],
    'subtotal_colectomy_with_end_ileostomy':[],
    'subtotal_colectomy_with_ileorectal_anastomosis':[],
    'intersphincteric_proctocolectomy':[],
    'intestinal_bypass':[],
    'proctectomy_with_ileal_pouch_anal_anastomosis':[],
    'proctosigmoidectomy_with_ileal_pouch_anal_anastomosis':[],
    'small_bowel_resection':[],
    'stricturoplasty':[],
    'anastomotic_resection_with_redo_anastomosis':[],
    'multiple_visceral_resection':[],
    'Appendectomy':[],
    'Cholecystectomy':[],
    'Cystoprostatectomy_Total_cystectomy':[],
    'Distal_sacrectomy':[],
    'Hysterectomy':[],
    'Inguinal_Ventral_hernia_repair':[],
    'Liver_resection':[],
    'Partial_cystectomy':[],
    'Partial_vaginectomy':[],
    'Plastics_flap_closure':[],
    'Salpingo_oophorectomy':[]}

    percentage = 0
    pt_cnt = 0
    test = []
    running_fxn(20,percentage)
    for pt in df.patient_id:
        pt_cnt +=1
        
        df_pt = df[df.patient_id==pt]

        if df_pt.shape[0]>1:
            pass
        else:

            pt_sx_list = df_pt.columns[df_pt.iloc[0,]==1].values

            condensed_pt_sx_list = []
            for item in pt_sx_list:
                if item != 'patient_id':
                    condensed_sx = sx_ddf.unique[sx_ddf.name==item].values[0]
                    condensed_pt_sx_list.append(condensed_sx)

            if condensed_pt_sx_list == []:
                test.append(pt)
            
            for item in output_dict:
                if item == 'patient_id':
                    output_dict[item].append(pt)
                elif item in condensed_pt_sx_list:
                    output_dict[item].append(1)
                else:
                    output_dict[item].append(0)
        if round(pt_cnt/num_of_pts*100) != percentage:
            percentage = round(pt_cnt/num_of_pts*100)
            if percentage in range(0,101,5):
                running_fxn(20,percentage)
    print('pt cnt: {}'.format(pt_cnt))

    df_out = pd.DataFrame(output_dict)


    df_final = pd.merge(df,df_out,how='inner',on='patient_id')

    print(test)

    print('Saving data to excel and pickle file...')

    writer = pd.ExcelWriter('S:\ERAS\crdb_output.xlsx')
    df_out.to_excel(writer,'sheet1')
    writer.close()

    pd.to_pickle(df_out,'S:\ERAS\crdb_output.pickle')

    

def pt_query():
    df = pd.read_pickle('S:\ERAS\crdb_output.pickle')

    tpc_list = ['total_proctocolectomy_with_end_ileostomy','total_proctocolectomy_with_ileal_pouch_anal_anastomosis','total_proctocolectomy_with_ileal_pouch_anal_anastomosis','total_proctocolectomy','total_abdominal_colectomy_with_end_ileostomy','total_abdominal_colectomy_with_ileorectal_anastamosis','total_abdominal_colectomy_with_ileorectal_anastamosis','total_abdominal_colectomy']

    pt_to_include = []
    for pt in range(df.shape[0]):
        pt_sx_list = df.columns[df.iloc[pt,]==1].values
        for item in pt_sx_list:
            if item in tpc_list:
                pt_to_include.append(df.patient_id.iloc[pt])

    # print('Number of patients in the following list of procedures {} was {}'.format(tpc_list,len(pt_to_include)))
    # print(pt_to_include)

    return pt_to_include

# load_and_pickle('S:\ERAS\CR_all.xlsx')
# extract_sx_data()
# create_sx_values()
tpc_pt_list = pt_query()

# print(tpc_pt_list)

def tpc_data(file,tpc_pt_list):
    df = pd.read_pickle('S:\ERAS\cr_datebase.pickle')
    # redcap = df.redcap_event_name.unique().tolist()
    df = df[df.patient_id.isin(tpc_pt_list)]
    df_pt = df[df.patient_id==30]
    
    redcap = df_pt.redcap_event_name.unique().tolist()
    cnt = 0
    print(df_pt.shape[1])
    for event in redcap:
        # print(df_pt[df_pt.redcap_event_name==event].shape)
        
        cnt+=df_pt[df_pt.redcap_event_name==event].dropna(axis=1).shape[1]
    print(cnt)


# tpc_data('S:\ERAS\CR_all.xlsx',tpc_pt_list)

def condense_rows():
    df = pd.read_pickle('S:\ERAS\cr_datebase.pickle') #whole database pickle
    df_dict = pd.read_pickle('S:\ERAS\colorecatal_dict.pickle') #loads dictionary pickle

    #creates a list of column headers for each type of event
    baseline_form_list = ['demographics','patient_info','patient_hx','family_hx']
    baseline_headers = df_dict.name[df_dict.form_name.isin(baseline_form_list)].tolist()

    preop_form_list = ['preop_evaluation_1','preop_evaluation_2','preop_evaluation_3']
    preop_headers = df_dict.name[df_dict.form_name.isin(preop_form_list)].tolist()

    neoadjuvant_form_list = ['neoadjuvant_treatment']
    neoadjuvant_headers = df_dict.name[df_dict.form_name.isin(neoadjuvant_form_list)].tolist()

    surgery_form_list = ['surgery_a','surgery_b','surgery_c','surgery_d']
    surgery_headers = df_dict.name[df_dict.form_name.isin(surgery_form_list)].tolist()

    pathology_form_list = ['pathology_a','pathology_b','pathology_c','pathology_d']
    pathology_headers = df_dict.name[df_dict.form_name.isin(pathology_form_list)].tolist()

    adjuvant_form_list = ['adjuvant_treatment']
    adjuvant_headers = df_dict.name[df_dict.form_name.isin(adjuvant_form_list)].tolist()

    postop_comp_form_list = ['post_op_complications_a','post_op_complications_b','post_op_complications_c','post_op_complications_d']
    postop_comp_headers = df_dict.name[df_dict.form_name.isin(postop_comp_form_list)].tolist()

    event_dict = {'baseline':baseline_headers,'preop':preop_headers,'neoadjuvant':neoadjuvant_headers,'surgery':surgery_headers,'pathology':pathology_headers,'adjuvant':adjuvant_headers,'postop_comp':postop_comp_headers}
    
    df_pt = df[df.patient_id==1] #just for testing. will replace with loop
    # regex_baseline = re.compile(r'baseline')
    # print(df_pt.redcap_event_name.tolist())
    # test = re.search(regex_baseline,df_pt.redcap_event_name.tolist())

    # test_df = df_pt[re.search(regex_baseline,df_pt.redcap_event_name)==1]

    # df_pt_baseline = df_pt[df_pt.redcap_event_name.str.contains('baseline')]
    # df_pt_baseline[]

    # for event in event_dict:
    #     print(df_pt[event_dict[event]].shape)


    # for event in event_dict:
    #     print(event_dict[event])




#need to convert data dictionary from the redcap database variable/field names to column names used in database output
#also need the form name which will be used for grouping the rows and condense the pts data to 1 row per pt

#will need to update the dictionary at some point
#there is an issue with _complete as it is not clear where these columns come from ? section headers. May just ignore for now




#baseline arms: demographics, patient info, family hx
#preop visit armS: preop eval 1,2,3
#neoadj arms: neoadjuvant_treatment
#surgery arms: surgery a,b,c,d
#pathology arms: pathology a,b,c,d
#adjuvant arms: adjuvant treatment
#post op comp amrs: post op complications a,b,c,d


def create_data_dict(input_dict):
    # df = pd.read_excel('S:\ERAS\sx_list_imput.xlsx')
    df = pd.read_excel(input_dict)
    main_output = [] #will be column name in CR database
    description_output = [] #description of column
    score = []
    unique_list = [] #unique names
    unique_code = [] #unique codes for each procedure
    form_name_list = []

    df_unique = pd.read_excel('S:\ERAS\CR_unique_dict.xlsx')
    procedure_dict = df_unique.to_dict()

    #iterate over all rows
    for row in df.iterrows():
        input_name = row[1].values[0] #gets main name
        form_name = row[1].values[3] #gest form name
        field_type = row[1].values[5] #gest field type
        text_names = row[1].values[7] #gets long string with values separated by "|", except for a few that have no string
        
        if field_type == 'checkbox':

            #for cells that have a strings separated by "|"
            try:
                text_list = text_names.split(' | ')
                for item in text_list:
                    match = False
                    for procedure in procedure_dict:      
                        find_procedure = re.search(procedure,item)
                        if find_procedure is not None:
                            unique_list.append(procedure_dict[procedure][0])
                            unique_code.append(procedure_dict[procedure][1])
                            score.append(procedure_dict[procedure][2])
                            form_name_list.append(form_name)
                            match = True
                            break #if match no need to look further
                    #checks to make sure there was a match
                    if not match:
                        unique_list.append('None')
                        unique_code.append(-1)
                        score.append(-1)
                        form_name_list.append(form_name)
                    regex = re.search(r'(\w+).*',item)
                    number = regex.group(1) #number value from string
                    procedure = regex.group(0) #whole string
                    main_output.append('{}___{}'.format(input_name,number)) #creates the unique name for each value in second column separted by "|"
                    description_output.append(procedure)
            
            #for cells that do not have a string in the second column
            except:
                item = input_name
                match = False
                for procedure in procedure_dict:      
                    find_procedure = re.search(procedure,item)
                    if find_procedure is not None:
                        unique_list.append(procedure_dict[procedure][0])
                        unique_code.append(procedure_dict[procedure][1])
                        score.append(procedure_dict[procedure][2])
                        form_name_list.append(form_name)
                        match = True
                        break #if match no neeed to look further
                #checks to make sure there was a match
                if not match:
                    unique_list.append('None')
                    unique_code.append(-1)
                    score.append(-1)
                    form_name_list.append(form_name)
                regex = re.search(r'(\w+).*',item)
                number = regex.group(1)
                procedure = regex.group(0)
                main_output.append(input_name)
                description_output.append(procedure)
        else: 
            item = input_name
            match = False
            for procedure in procedure_dict:      
                find_procedure = re.search(procedure,item)
                if find_procedure is not None:
                    unique_list.append(procedure_dict[procedure][0])
                    unique_code.append(procedure_dict[procedure][1])
                    score.append(procedure_dict[procedure][2])
                    form_name_list.append(form_name)
                    match = True
                    break #if match no neeed to look further
            #checks to make sure there was a match
            if not match:
                unique_list.append('None')
                unique_code.append(-1)
                score.append(-1)
                form_name_list.append(form_name)
            regex = re.search(r'(\w+).*',item)
            number = regex.group(1)
            procedure = regex.group(0)
            main_output.append(input_name)
            description_output.append(procedure)
        
    df_out = pd.DataFrame(main_output,columns=['name'])
    df_out['score'] = score
    df_out['description'] = description_output
    df_out['unique'] = unique_list
    df_out['code'] = unique_code
    df_out['form_name'] = form_name_list
    writer = pd.ExcelWriter('S:\ERAS\colorectal_dict.xlsx')
    df_out.to_excel(writer,'Sheet1')
    writer.close()

    pd.to_pickle(df_out,'S:\ERAS\colorecatal_dict.pickle')

# create_data_dict('S:\ERAS\colorectal_registry_dictionary_06012017.xlsx')


def testing():
    df = pd.read_pickle('S:\ERAS\cr_datebase.pickle')
    df_dict = pd.read_pickle('S:\ERAS\sx_list_dict_comp_test.pickle')
    df_col_list = df.columns.tolist()
    df_dict_list = df_dict.name.tolist()
    print(list(set(df_col_list)-set(df_dict_list)))
    # print(list(set(df_dict_list)-set(df_col_list)))

# testing()
condense_rows()

"""
1 - no
2- proctosigmoidectomy
3,4,5,6, - no
7,8,9 - other
9 - proctosigmoidectomy



"""