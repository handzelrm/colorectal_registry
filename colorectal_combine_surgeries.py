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

    # sx_no_dict = ['prim_sx_colonca_a___35', 'prim_sx_colonca_a___36', 'prim_sx_colonpolyp_a___34', 'prim_sx_colonpolyp_a___35', 'sx_anastomosis_a', 'sx_anastamosis_ibd_a', 'sx_comb_service_a___16', 'sx_comb_service_a___17', 'sx_comb_service_a___18', 'sx_comb_service_a___19', 'sx_comb_service_a___20', 'sx_comb_service_a___21']

    # for item in sx_list:
    #     # print(item)
    #     print(item)
    #     if item == 'sx_comb_service_a___17':
    #         print(item)
    #     if item in sx_no_dict:
    #         sx_list.remove(item)
    #         # print(item)
    #         # print('match')
    # # print(sx_list)
    # return
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
        

    # sx_no_dict = ['prim_sx_colonca_a___35', 'prim_sx_colonca_a___36', 'prim_sx_colonpolyp_a___34', 'prim_sx_colonpolyp_a___35', 'sx_anastomosis_a', 'sx_anastamosis_ibd_a', 'sx_comb_service_a___16', 'sx_comb_service_a___17', 'sx_comb_service_a___18', 'sx_comb_service_a___19', 'sx_comb_service_a___20', 'sx_comb_service_a___21']


    # return

    # print(df.head())
    # return
    # print(test)
    # print(len(test))
    # return

    # print(df.columns)

    percentage = 0
    pt_cnt = 0
    test = []
    running_fxn(20,percentage)
    for pt in df.patient_id:
        # cnt = 0
        pt_cnt +=1
        # print(pt_cnt)
        
        df_pt = df[df.patient_id==pt]

        #4.4 seconds!
        if df_pt.shape[0]>1:
            pass
        else:
            # print(df_pt.shape)
            pt_sx_list = df_pt.columns[df_pt.iloc[0,]==1].values
            # output_dict['patient_id'].append(pt)
            # if 'patient_id' in pt_sx_list:
            #     del pt_sx_list[0] #del patient_id from list

            # print(pt_sx_list)
            # print(pt_sx_list)
            # print(pt_sx_list)
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
                    # print('match')
                else:
                    output_dict[item].append(0)
        # return
        """
        #851.1s
        for col in df_pt.columns:
            if df_pt.shape[0]>1:
                pass
            else:
                if col == 'patient_id':
                    pass
                else:
                    try:
                        sx_column = sx_ddf.unique[sx_ddf.name==col].values[0]
                    except:
                        if col not in test:
                            test.append(col)
                        
                        
                    if df_pt[col].values[0] == 1:
                        # cnt +=1                        
                        output_dict[sx_column].append(1)
                    else:
                        output_dict[sx_column].append(0)
                        """
                        
                        # return col
                    # print(col)
                    # print(df_pt[col].values[0])
        # cnt_list.append(cnt)

        #display percentage
        if round(pt_cnt/num_of_pts*100) != percentage:
            percentage = round(pt_cnt/num_of_pts*100)
            if percentage in range(0,101,5):
                running_fxn(20,percentage)
    # print(test)
    print('pt cnt: {}'.format(pt_cnt))
    # print(max(cnt_list))
    # for sx in output_dict:
    #     print(len(output_dict[sx]))
    df_out = pd.DataFrame(output_dict)
    # print(df_out.patient_id)

    df_final = pd.merge(df,df_out,how='inner',on='patient_id')

    print(test)

    # print(df.shape)
    # print(df_out.shape)
    # print(df_final.shape)

    print('Saving data to excel and pickle file...')

    writer = pd.ExcelWriter('S:\ERAS\crdb_output.xlsx')
    df_out.to_excel(writer,'sheet1')
    writer.close()

    pd.to_pickle(df_out,'S:\ERAS\crdb_output.pickle')

    

def pt_query():
    df = pd.read_pickle('S:\ERAS\crdb_output.pickle')

    tpc_list = ['total_proctocolectomy_with_end_ileostomy','total_proctocolectomy_with_ileal_pouch_anal_anastomosis','total_proctocolectomy_with_ileal_pouch_anal_anastomosis','total_proctocolectomy','total_abdominal_colectomy_with_end_ileostomy','total_abdominal_colectomy_with_ileorectal_anastamosis','total_abdominal_colectomy_with_ileorectal_anastamosis','total_abdominal_colectomy']

    # # print(df.patient_id[df==1])
    # print(df.shape[0])
    # print(df.columns[df.iloc[0,]==1].values)
    # print(df.head())
    # print(df.head())
    pt_to_include = []
    for pt in range(df.shape[0]):
        pt_sx_list = df.columns[df.iloc[pt,]==1].values
        for item in pt_sx_list:
            if item in tpc_list:
                pt_to_include.append(df.patient_id.iloc[pt])

    print('Number of patients in the following list of procedures {} was {}'.format(tpc_list,len(pt_to_include)))
    print(pt_to_include)







load_and_pickle('S:\ERAS\CR_all.xlsx')
extract_sx_data()
create_sx_values()
pt_query()

"""
1 - no
2- proctosigmoidectomy
3,4,5,6, - no
7,8,9 - other
9 - proctosigmoidectomy



"""