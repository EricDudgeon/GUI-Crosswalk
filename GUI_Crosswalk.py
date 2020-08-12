# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 17:56:47 2020

@author: Eric Dudgeon
"""
import textwrap
from tkinter import *
#from tkinter import ttk
from tkinter import filedialog
import pandas as pd
import numpy as np
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import re
import sys
import os
from Levenshtein import *


csv_crosswalk = "Crosswalk Template.csv"
csv_fomulary = "MHSGFormulary.csv"
dosage_dict = pd.read_csv("Dictionary List.csv")

def run_pyx():
    cw = pd.read_csv(csv_crosswalk, encoding="utf-8")
    cw.index.name = None
    extract = pd.read_csv(csv_fomulary, low_memory=False, encoding="utf-8")
    extract.index.name = None
    extract_org = pd.DataFrame(extract, columns = ["LABEL_DESC","ADM_ID","GENERIC_NAME","BRAND_NAME","FORM","DOSE","VOLUME"])
    
    # #crosswalk columns upper
    cw["Client Med ID"] = cw["Med_ID"]
    cw["LABEL_DESC"] = cw["Label_Desc"].str.upper().astype("string")
    cw["BRAND_NAME"] = cw["Brand_Name"].str.upper().astype("string")
    dosage_form_original = cw["Dosage_Form"].astype("string")
    cw["FORM"] = cw["Dosage_Form"].str.upper().astype("string")
    cw["DOSE"]= cw["Strength"].str.upper().astype("string")
    #cw["GENERIC_NAME"] = cw["PureGenericName"].str.upper().astype("string")
    cw["VOLUME"] =cw["Volume"].str.upper().astype("string")
    
    #Extract columns upper
    extract["LABEL_DESC"] = extract["LABEL_DESC"].str.upper().astype("string")
    extract["BRAND_NAME"] = extract["BRAND_NAME"].str.upper().astype("string")
    extract["FORM"] = extract["FORM"].str.upper().astype("string")
    extract["DOSE"] = extract["DOSE"].str.upper().astype("string")
    extract["GENERIC_NAME"] = extract["GENERIC_NAME"].str.upper().astype("string")
    
    #modifying data to only include useful rows, and defining variables to use
    cw = pd.DataFrame(cw, columns=["LABEL_DESC","BRAND_NAME","FORM","DOSE","VOLUME", "Client Med ID"])
    extract = pd.DataFrame(extract, columns = ["LABEL_DESC","ADM_ID","GENERIC_NAME","BRAND_NAME","FORM","DOSE","VOLUME"])
    
    
    # Creating a dictionary of dosage form corrections
    # Dictionary created in seperate module by using the unique values method on both the formulary extact and crosswalk formularies
    # from all of wave travis and pendleton (extract formulary has 108 and cross walk has 230 unique dosage forms)
    # abbrev and errors linked to correct Formulary extract dosage forms
    list1 = dosage_dict.values.tolist()
    dictionary = {}
    for each in list1:
        dictionary[each[0]] = each[1]
    
    # Compares cross formulary against the dictionary and replaces any errors
    # dosage_form_corrected = cw['FORM'].str.upper().replace("-"," ", regex = True).map(dictionary)
    # cw["FORM"] = dosage_form_corrected
    dosage_form_corrected = cw["FORM"]
    forms_clean = []
    for form in dosage_form_corrected:
        form = form.replace("-"," ")
        if form in dictionary:
            forms_clean.append(dictionary[form])
            dosage_form_corrected = forms_clean
        else:
            forms_clean.append(form)
            dosage_form_corrected = forms_clean
    
    dosage_form_corrected = pd.Series(dosage_form_corrected)
    
    extract1 = pd.DataFrame(extract)
    extract1.drop_duplicates(inplace = True)
    
    replacement_dict = {" CR ":" CR RELEASE ", " ER ":" EXTENDED RELEASE ", "1 G": "1000 MG", "U/D": "UD", " NS": " SODIUM CHLORIDE", 
                         "D5W":"DEXTROSE 5% WATER","VIAL":" VIAL INJECTION ","-":" "," 1G ":"1000 MG","2 G": "2000 MG", " EC ":" ENTERIC COATED ",
                       " INJ ":" INJECTION ", "SYRINGE":" SYRINGE INJECTION "," TAB-CHEW ":"TABLET CHEWABLE"," D5 ":"Dextrose 5%","*LOOK-ALIKE/SOUND-ALIKE*":"","**HIGH ALERT**":"", " **HA**":"","*":"",
                       "*HAZARDOUS MEDICATION*":"","*SINGLE DOSE*":""}
    labels_conversion1 = {"TABLET":" TAB ", "CAPSULE": " CAP ","TAB,":" TAB ", "CAP,":" CAP ", "INJECTION": " INJ ",
                          "SUSPENSION":" SUSP ","SOLUTION": " SOL "," TAB":" TAB ", " CAP":" CAP "," INJ":" INJ "," SOL":" SOL "," SUSP":" SUSP "}
    labels_conversion2 = {" TAB ":" TABLET "," CAP ": " CAPSULE "," INJ ":" INJECTION "," SUSP ":" SUSPENSION "," SOL ":" SOLUTION " }
    
    
    ###Crosswalk concat 3
    concat3 = cw['LABEL_DESC'].str.replace(" ","",).map(str) + cw['DOSE'].fillna("").str.replace(" ","").map(str) + dosage_form_corrected.str.replace(" ","")
    concat3 = concat3.astype("string")
    concat3 = concat3.apply(str)
    
    ###crosswalk concat 1
    concat1 = (cw['LABEL_DESC'].map(str) + " " + cw['BRAND_NAME'].fillna("").map(str) + " " + cw['DOSE'].fillna("").str.replace(" ","").map(str) + " " + 
               dosage_form_corrected + " " + cw["VOLUME"].fillna("").str.replace(" ","").map(str)) 
    concat1 = concat1.astype("string")
    concat1 = concat1.apply(str)
    
    ### Label Desc Conv
    label_desc = extract1['LABEL_DESC'].astype("string")
    
    
    def multipleReplace(text, wordDict):
        for key in wordDict:
            text = text.replace(key, wordDict[key])
        return text
    
    concat1 = concat1.apply(lambda x: multipleReplace(x,replacement_dict))
    label_desc = label_desc.apply(lambda x: multipleReplace(x,replacement_dict))
    
    
    
    label_desc = label_desc.apply(lambda x: multipleReplace(x,labels_conversion1))
    label_desc = label_desc.apply(lambda x: multipleReplace(x,labels_conversion2))
    concat1 = concat1.apply(lambda x: multipleReplace(x,labels_conversion1))
    concat1 = concat1.apply(lambda x: multipleReplace(x,labels_conversion2))
    label_desc = label_desc.str.replace(' MG', "MG").str.replace("  ", " ").str.replace("   "," ").str.replace(" ML","ML")
    
    
    def lev_equation_two(list_a, list_b, item_list, conf_list, scorer):
        for item in list_a:
            confidence = process.extractOne(item, list_b, scorer = scorer)
            item_list.append(item)
            conf_list.append(confidence)
    
    item_list = []
    confidence_list = []
    lev_equation_two(concat3, label_desc, item_list, confidence_list, fuzz.ratio)
    
    
    match_list = []
    score_list = []
    index_match = []
    
    for match, score, row in confidence_list:
        match_list.append(match)
        score_list.append(score)
        index_match.append(row)
    
    
    index = pd.Series(index_match)
    results_list = []
    for item in index:
        results = extract_org.iloc[item]
        results_list.append(results)
        
    label_final = []
    adm_final = []
    for label, adm, generic, brand, form, dose, volume in results_list:
        label_final.append(label)
        adm_final.append(adm)
    
    
    
    item_list2 = []
    confidence_list2 = []
    match_list2 = []
    score_list2 = []
    index_match2 = []
    lev_equation_two(concat1, label_desc, item_list2, confidence_list2, fuzz.token_set_ratio)
    for match, score, row in confidence_list2:
        match_list2.append(match)
        score_list2.append(score)
        index_match2.append(row)
        
    index2 = pd.Series(index_match2)
    results_list2 = []
    for item in index2:
        results = extract_org.iloc[item]
        results_list2.append(results)
        
    label_final2 = []
    adm_final2 = []
    for label, adm, generic, brand, form, dose, volume in results_list2:
        label_final2.append(label)
        adm_final2.append(adm)
    
    
    ##creates output dataframe
    output = pd.DataFrame(cw, columns = ['LABEL_DESC',"BRAND_NAME","Client Med ID"])
    output.rename( columns={ "LABEL_DESC" :'Client Generic Name', "BRAND_NAME":"Client Brand Name"}, inplace=True )
    output["Client Brand Name"] = output["Client Brand Name"].fillna("Not Given")
    output["Client Strength"] = cw["DOSE"].fillna("Not Given")
    output["Client Dosage Form"] = pd.Series(dosage_form_original)
    output["Client Volume"] = cw["VOLUME"].fillna("Not Given")
    output["Cerner Label Desc Proposed 1"] = pd.Series(label_final)
    output["ADM Proposed 1"] = pd.Series(adm_final)
    output["Confidence Rating 1"] = pd.Series(score_list)
    output["Cerner Label Desc Proposed 2"] = pd.Series(label_final2)
    output["ADM Proposed 2"] = pd.Series(adm_final2)
    output["Confidence Rating 2"] = pd.Series(score_list2)
    output["Label Final"] = ""
    output["ADM Final"] = ""
    output["Confidence Final"] = ""
    output["Concat1"] = pd.Series(concat1)
    
    
    def func(row):
        if row['ADM Proposed 1'] == row['ADM Proposed 2']:
            return row["Cerner Label Desc Proposed 1"]
            
        elif row['Confidence Rating 1'] > 90 and "SYNTHROID" not in row["Client Brand Name"]:
            return row["Cerner Label Desc Proposed 1"]
            
        else:
            return row["Cerner Label Desc Proposed 2"]
    
    
    def func1(row):
        if row['ADM Proposed 1'] == row['ADM Proposed 2']:
            return row["ADM Proposed 1"]
            
        elif row['Confidence Rating 1'] > 90 and "SYNTHROID" not in row["Client Brand Name"]:
            return row["ADM Proposed 1"]
            
        else:
            return row["ADM Proposed 2"]
    
    
    def func2(row):
        if row['ADM Proposed 1'] == row['ADM Proposed 2']:
            return (row["Confidence Rating 1"] + row["Confidence Rating 2"]) /2
            
        elif row['Confidence Rating 1'] > 90 and "SYNTHROID" not in row["Client Brand Name"]:
            return row["Confidence Rating 1"]
            
        else:
            return (row["Confidence Rating 1"] + row["Confidence Rating 2"]) /2
    
    
    output['Label Final'] = output.apply(func, axis=1)
    output['ADM Final'] = output.apply(func1, axis=1)
    output["Confidence Final"] = output.apply(func2, axis=1)
    
    
    output_final = pd.DataFrame(output, columns=["Client Generic Name","Client Brand Name","Client Med ID", "Client Strength","Client Dosage Form","Client Volume","Label Final", "ADM Final", "Confidence Final",])
    output_final.rename( columns={ "Label Final" :'Cerner Label Desc Proposed', "ADM Final":"ADM Proposed", "Confidence Final":"Confidence"}, inplace=True )
    output_final["Audit"] = ""
    
    
    l = output_final["Cerner Label Desc Proposed"].str.lower().str.replace(",","").str.replace("/"," ").str.replace("-"," ").str.replace("["," ").str.replace("]"," ").str.split()
    s = output_final["Client Strength"].str.lower().str.replace(",","").str.split()
    v = output_final["Client Volume"].str.lower().str.split()
    n = output_final["Client Generic Name"].str.replace("/"," ").str.replace("-"," ").str.replace("["," ").str.replace("]"," ").str.split()
    volume_audit = []
    strength_audit = []
    name_audit = []
    for each in n:
        name = each[0]
        name_audit.append(name)
    for each in v:
        volume = each[0]
        volume_audit.append(volume)
    for each in s:
        strength = each[0]
        strength_audit.append(strength) 
    ls = pd.Series(l)
    df = pd.DataFrame(ls)
    df["Name"] = pd.Series(name_audit)
    df["Name"] = df["Name"].str.lower()
    df["Strength"] = pd.Series(strength_audit)
    df["Volume"] = pd.Series(volume_audit)
    df["Name Int"] = ""
    df["Modified Release"] =""
    df["Form"] = output_final["Client Dosage Form"].str.replace("-"," ")
    df["Brand"] = output_final["Client Brand Name"]
    df["Full Name"] = output_final["Client Generic Name"].str.replace("/"," ").str.replace("-"," ").str.replace("["," ").str.replace("]"," ")
    df["Audit"] = ""
    
    
    df["Label"] = output_final["Cerner Label Desc Proposed"]
    df["Name Int"] = df["Full Name"].str.extract('(\d+%)').fillna("no")
    df["proposed Int"] = df["Label"].str.extract('(\d+%)').fillna("no")
    
    
    def match_func(row):
        a = row["Name Int"]
        b = row["proposed Int"]
        if b == a:
            return "yes"
        else:
            return "no"
    
    
    df["Match"] = df.apply(match_func, axis=1)
    
    
    
    def release_audit(row):
        if " ER" in row["Form"] or " ER" in row["Brand"] or "ER " in row["Form"]:
            return "ER"
        elif " DR" in row["Form"] or " DR" in row["Brand"] or "DR " in row["Form"]:
            return "DR"
        elif " CR" in row["Form"] or " CR" in row["Brand"] or "CR " in row["Form"]:
            return "CR"
        elif " XL" in row["Form"] or " XL" in row["Brand"] or "XL " in row["Form"]:
            return "XL"
        elif " CR" in row["Form"] or " CR" in row["Brand"] or "CR " in row["Form"]:
            return "CR"
        elif " XR" in row["Form"] or " XR" in row["Brand"] or "XR " in row["Form"]:
            return "XR"
        elif " SR" in row["Form"] or " SR" in row["Brand"] or "SR " in row["Form"]:
            return "SR"
        elif " 24 hr" in row["Form"] or " 24 hr" in row["Brand"]:
            return "24 hr"
        elif " 12 hr" in row["Form"] or " 12 hr" in row["Brand"]:
            return "12 hr"
        elif " EC" in row["Form"] or " EC" in row["Brand"] or "EC " in row["Form"]:
            return "EC"
        elif " Chew" in row["Form"] or " Chew" in row["Brand"] or"Chew" in row["Form"] or "Chewable" in row["Form"]:
            return "Chew"
        elif "Disintegrating" in row["Form"] or "Disintegrating" in row["Brand"]:
            return "Disintegrating"
        else:
            return "no"
    
    
    df["Modified Release"] = df.apply(release_audit,axis=1)
    
    
    def audit_func(row):
        l = row["Cerner Label Desc Proposed"]
        v = row["Volume"]
        s = row["Strength"]
        mr = row["Modified Release"]
        n = row["Name"]
        fn = row["Full Name"]
        z = row["Match"]
        if v in l and s in l and n in l and mr == "no" and z == "yes":
            return "OK"
        elif v in l and s == "not" and n in l and mr =="no" and z == "yes":
            return "OK"
        elif v == "not" and s in l and n in l and mr =="no" and z== "yes":
            return "OK"
        elif v in l and s in l and n in l and mr != "no" and z=="yes":
            return "OK if "+ mr +" or equiv"
        elif v in l and s == "not" and n in l and mr !="no" and z== "yes":
            return "OK if "+ mr +" or equiv"
        elif v == "not" and s in l and n in l and mr !="no" and z=="yes":
            return "OK if "+ mr +" or equiv"
        elif v in l and s not in l and n in l and mr =="no" and z=="yes":
            return "Check Strength"
        elif v not in l and s in l and n in l and mr=="no" and z=="yes":
            return "Check Volume"
        elif v in l and s in l and n not in l and mr=="no" and z=="yes":
            return "Check Med"
        elif v in l and s not in l and n in l and mr !="no" and z=="yes":
            return "Check Strength and " + mr
        elif v not in l and s in l and n in l and mr!="no" and z=="yes":
            return "Check Volume and " + mr
        elif v in l and s in l and n not in l and mr!="no" and z=="yes":
            return "Check Med and " + mr
        elif v == "not" and s == "not":
            return "Unknown"
        else:
            return "Review"
    
    df["Audit"] = df.apply(audit_func, axis=1)
    
    
    output_final["Audit"] = df["Audit"]
    
    
    def lev_func_three(list_a, list_b, item_list, conf_list):
        for item in list_a:
            confidence = process.extract(item, list_b, limit=3, scorer = fuzz.token_set_ratio)
            item_list.append(item)
            conf_list.append(confidence)
    
    
    item_list3 = []
    confidence_list3 = []
    match_list3 = []
    score_list3 = []
    index_match3 = []
    lev_func_three(concat1, label_desc, item_list3, confidence_list3)
    
    
    options_list = []
    item_list_final = []
    brand_options = []
    brand = cw["BRAND_NAME"].fillna("Not Given")
    options_names = cw["LABEL_DESC"].map(str) + " " + cw["DOSE"].fillna("").map(str) + " " + dosage_form_original.fillna("") + " " + cw["VOLUME"].fillna("").map(str)
    medication_ID = []
    med = output["Client Med ID"].fillna("Not Given")
    for item in med:
        medication_ID.append(item)
        medication_ID.append("---")
        medication_ID.append("---")
        
    for item in brand:
        brand_options.append(item)
        brand_options.append("---")
        brand_options.append("---")
    
    for item in options_names:
        item_list_final.append(item)
        item_list_final.append("---")
        item_list_final.append("---")
    for r1, r2, r3 in confidence_list3:
        options_list.append(r1)
        options_list.append(r2)
        options_list.append(r3)
        
    for match, score, row in options_list:
        match_list3.append(match)
        score_list3.append(score)
        index_match3.append(row)
        
    index3 = pd.Series(index_match3)
    results_list3 = []
    for item in index3:
        results = extract_org.iloc[item]
        results_list3.append(results)
        
    label_final3 = []
    adm_final3 = []
    for label, adm, generic, brand, form, dose, volume in results_list3:
        label_final3.append(label)
        adm_final3.append(adm)
    
    
    options_list_output = pd.DataFrame()
    options_list_output["Client Generic Name"] = pd.Series(item_list_final)
    options_list_output["Client Brand Name"] = pd.Series(brand_options)
    options_list_output["Client Med ID"] = pd.Series(medication_ID)
    options_list_output["Addl Label Options"] = pd.Series(label_final3)
    options_list_output["ADM_ID"] = pd.Series(adm_final3)
    
    with pd.ExcelWriter("Crosswalk Automation Output.xlsx") as writer:
        output_final.to_excel(writer,sheet_name = "Crosswalk Template")
        options_list_output.to_excel(writer,sheet_name = "Options", index=None)
    comp_label.grid(row=4, pady=5,padx=20, sticky=W)
##################################   
    
     
def run_labels():
    
    cw = pd.read_csv('Crosswalk Template.csv', encoding="utf-8")
    cw.index.name = None
    extract = pd.read_csv('MHSGFormulary.csv', low_memory=False, encoding="utf-8")
    extract.index.name = None
    extract_org = pd.DataFrame(extract, columns = ["LABEL_DESC","ADM_ID","GENERIC_NAME","BRAND_NAME","FORM","DOSE","VOLUME"])
    
    
    list1 = dosage_dict.values.tolist()
    dictionary = {}
    for each in list1:
        dictionary[each[0]] = each[1]
    
    replacement_dict = {" CR ":" CR RELEASE ", " ER ":" EXTENDED RELEASE ", "1 G": "1000 MG", "U/D": "UD", " NS": " SODIUM CHLORIDE", 
                         "D5W":"DEXTROSE 5% WATER","VIAL":" VIAL INJECTION ","-":" "," 1G ":"1000 MG","2 G": "2000 MG", " EC ":" ENTERIC COATED ",
                       " INJ ":" INJECTION ", "SYRINGE":" SYRINGE INJECTION "," TAB-CHEW ":"TABLET CHEWABLE"," D5 ":"Dextrose 5%","*LOOK-ALIKE/SOUND-ALIKE*":"","**HIGH ALERT**":"", " **HA**":"","*":"",
                       "/1ML":"/ML","1GM":"1000 MG", "2GM":"2000 MG","*SINGLE DOSE*":"","*HAZARDOUS MEDICATION*":""}
    
    labels_conversion1 = {"TABLET":" TAB ", "CAPSULE": " CAP ","TAB,":" TAB ", "CAP,":" CAP ", "INJECTION": " INJ ",
                          "SUSPENSION":" SUSP ","SOLUTION": " SOL "," TAB":" TAB ", " CAP":" CAP "," INJ":" INJ "," SOL":" SOL "," SUSP":" SUSP "}
    labels_conversion2 = {" TAB ":" TABLET "," CAP ": " CAPSULE "," INJ ":" INJECTION "," SUSP ":" SUSPENSION "," SOL ":" SOLUTION " }
    
    
    
    # #crosswalk columns upper
    
    cw["OmniDesc"] =cw["Label_Desc"].str.upper().astype("string")
    
    #Extract columns upper
    extract["LABEL_DESC"] = extract["LABEL_DESC"].str.upper().astype("string")
    extract["BRAND_NAME"] = extract["BRAND_NAME"].str.upper().astype("string")
    extract["FORM"] = extract["FORM"].str.upper().astype("string")
    extract["DOSE"] = extract["DOSE"].str.upper().astype("string")
    extract["GENERIC_NAME"] = extract["GENERIC_NAME"].str.upper().astype("string")
    
    #modifying data to only include useful rows, and defining variables to use
    cw = pd.DataFrame(cw, columns=["OmniDesc"])
    extract = pd.DataFrame(extract, columns = ["LABEL_DESC","ADM_ID","GENERIC_NAME","BRAND_NAME","FORM","DOSE","VOLUME"])
    
    
    
    extract1 = pd.DataFrame(extract)
    extract1.drop_duplicates(inplace = True)
    
    
    ###Crosswalk concat 3
    concat3 = cw['OmniDesc'].map(str)
    concat3 = concat3.astype("string")
    concat3 = concat3.apply(str)
    
    ###crosswalk concat 1
    concat1 = cw['OmniDesc'].map(str)
    concat1 = concat1.astype("string")
    concat1 = concat1.apply(str)
    
    ### Label Desc Conv
    label_desc = extract1['LABEL_DESC'].astype("string")
    
    
    def multipleReplace(text, wordDict):
        for key in wordDict:
            text = text.replace(key, wordDict[key])
        return text
    
    concat3 = concat3.apply(lambda x: multipleReplace(x,replacement_dict))
    label_desc = label_desc.apply(lambda x: multipleReplace(x,replacement_dict))
    
    
    
    label_desc = label_desc.apply(lambda x: multipleReplace(x,labels_conversion1))
    label_desc = label_desc.apply(lambda x: multipleReplace(x,labels_conversion2))
    concat3 = concat3.apply(lambda x: multipleReplace(x,labels_conversion1))
    concat3 = concat3.apply(lambda x: multipleReplace(x,labels_conversion2))
    # label_desc = label_desc.apply(lambda x: multipleReplace(x,dictionary))
    # concat3 = concat3.apply(lambda x: multipleReplace(x,dictionary))
    label_desc = label_desc.str.replace(" ", "")
    concat3 = concat3.str.replace(" ","")
    
    
    def lev_equation_two(list_a, list_b, item_list, conf_list, scorer):
        for item in list_a:
            confidence = process.extractOne(item, list_b, scorer = scorer)
            item_list.append(item)
            conf_list.append(confidence)
    
    
    item_list = []
    confidence_list = []
    lev_equation_two(concat3, label_desc, item_list, confidence_list, fuzz.ratio)
    
    
    match_list = []
    score_list = []
    index_match = []
    
    for match, score, row in confidence_list:
        match_list.append(match)
        score_list.append(score)
        index_match.append(row)
    
    
    index = pd.Series(index_match)
    results_list = []
    for item in index:
        results = extract_org.iloc[item]
        results_list.append(results)
        
    label_final = []
    adm_final = []
    for label, adm, generic, brand, form, dose, volume in results_list:
        label_final.append(label)
        adm_final.append(adm)
    
    
    
    output = pd.DataFrame(cw, columns = ["OmniDesc"])
    output["Cerner Label Desc Proposed"] = pd.Series(label_final)
    output["ADM Proposed"] = pd.Series(adm_final)
    output["Confidence Rating"] = pd.Series(score_list)
    
    
    f = open("Crosswalk Automation Output.xlsx","w")
    output.to_excel("Crosswalk Automation Output.xlsx")
    comp_label.grid(row=4, pady=5,padx=20, sticky=W)
    
    

def information():
    info_window = Toplevel()
    info_window.title("Help Menu")
    #info_window.iconbitmap("icon.ico")
    info_frame = Frame(info_window)
    label = Label(info_window,text='"Run Only Label Desc" required feilds in row one of template:')
    label.grid(row=0, pady=5,padx=10, sticky=W)
    label2 = Label(info_window, text="Label_Desc")
    label2.grid(row=1, pady=5,padx=10, sticky=W)
    label3 = Label(info_window, text='"Run Full Data" required feilds in row one of template:')
    label3.grid(row=2, pady=5,padx=10, sticky=W)
    label3 = Label(info_window, text="Label_Desc, Brand_Name, Med_ID, Dosage_Form, Strength, Volume")
    label3.grid(row=3, pady=5,padx=10, sticky=W)
    label4 = Label(info_window, text="Must use template provided in folder. Do not change name or file type until you save a copy")
    label4.grid(row=4, pady=5,padx=10, sticky=W)
    label5 = Label(info_window, text="Estimated Time to completion is 10 minutes per 2000 rows")
    label5.grid(row=5, pady=5,padx=10, sticky=W)

# def highlighter():
#     output = pd.read_excel('Crosswalk Automation Output.xlsx', encoding="utf-8")
#     output_final = pd.DataFrame(output)
#     output_final = output_final.iloc[:,1:11]
#     options_list = pd.read_excel('Crosswalk Automation Output.xlsx', encoding="utf-8", sheet_name="Options")
#     options_list_output = pd.DataFrame(options_list)
#     options_list_output = options_list_output.iloc[:1:6]
#     def highlight_output(row):        
#         n = row["Client Generic Name"]
#         if ("AMIODARONE" in n or "CARBAMAZEPINE" in n or "LEVOTHYROXINE" in n or "CLOZAPINE" in n or "CYCLOSPORINE" in n
#            or "DIVALPROEX SODIUM" in n or "DOFETILIDE" in n or "ETHOSUXIMIDE" in n or "LAMOTRIGINE" in n or "LEVETIRACETAM" in n
#            or "MESALAMINE" in n or "MYCOPHENOLATE" in n or "PHENYTOIN" in n or "TACROLIMUS" in n or "THEOPHYLLINE" in n or "THYROID DESICCATED" in n
#            or "TOPIRAMATE" in n or "WARFARIN" in n or "ZONISAMIDE" in n or "BUPIVACAINE" in n or "ROPIVACAINE" in n or "ADENOSINE" in n 
#            or "ARIPIPRAZOLE" in n or "ATROPINE" in n or "BELIMUMAB" in n or "DIPHENHYDRAMINE" in n or "EVOLOCUMAB" in n or "EXENATIDE" in n 
#            or "FENTANYL" in n or "GOLIMUMAB" in n or "IXEKIZUMAB" in n or "KETOROLAC" in n or "LORAZEPAM" in n or "MEDROXYPROGESTERONE" in n 
#            or "OCTREOTIDE" in n or "RANIBIZUMAB" in n):
#             return ['background-color: yellow'] * 10
#         else:
#             return ['background-color: none'] * 10
        
#     output_final = output_final.style.apply(highlight_output, axis=1)
    
#     with pd.ExcelWriter("Crosswalk Automation Output.xlsx") as writer:
#         output_final.to_excel(writer,sheet_name = "Crosswalk Template")
#         options_list_output.to_excel(writer,sheet_name = "Options", index=None)
    
####Main Window
root = Tk()
root.title("Crosswalk Automation")
root.geometry('{}x{}'.format(205, 200))
#root.iconbitmap("icon.ico")

# create all of the main containers
top_frame = Frame(root)
center = Frame(root)
# btm_frame = Frame(root)
# btm_frame2 = Frame(root, pady=3)

# layout all of the main containers
root.grid_rowconfigure(1, weight=1)
root.grid_columnconfigure(0, weight=1)


top_frame.grid(row=0,sticky=N+E+S+W)
center.grid(row=1,sticky=W+N+E+S, pady=10)
# center.grid(row=3,sticky=W+N+E+S,pady=10)
# btm_frame.grid(row=2, sticky=W+N+E+S, pady=10)
# btm_frame2.grid(row=4)


### CREATING BUTTONS
just_desc = Button(top_frame, text="Run Only Label Descriptions", command=run_labels)
pyxis = Button(top_frame, text="Run Full Data", command = run_pyx)
#button3 = Button(top_frame, text = "Highlight Narrow Therp", command = highlighter)

just_desc.grid(row=0,sticky=W+E+N+S, pady=5, padx=20)
pyxis.grid(row=1,sticky=W+E+S+N, pady=5, padx=20)
#button3.grid(row=2,sticky=W+E+S+N,pady=5,padx=20)
######CREATING RUNNING AND COMPLETE LABELS
run_label = Label(top_frame,text='Running...ETA 15 minutes')
#run_label.grid(row=2, pady=5,padx=20, sticky=W)
comp_label = Label(top_frame, text="Complete! Check Output file")
#comp_label.grid(row=3, pady=5,padx=20, sticky=W)

####CREATING MENUBAR
menu = Menu(root)
root.config(menu=menu)

subMenu = Menu(menu)
menu.add_cascade(label = "File", menu=subMenu)
subMenu.add_command(label="Exit", command=root.destroy)

updateMenu = Menu(menu)
menu.add_cascade(label = "Read me", menu=updateMenu)
updateMenu.add_command(label="Information", command=information)


root.mainloop()