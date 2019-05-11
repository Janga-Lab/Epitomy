from flask import Markup
from flask import Flask, request, render_template, jsonify, abort, make_response
from elasticsearch import Elasticsearch
from bs4 import BeautifulSoup
from human_genes import genes_coordinates
from mouse_genes import genes_coordinates_mouse
from multiprocessing.pool import ThreadPool
import requests, smtplib, urllib3, json
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def human_annotations(ensg,genename):
    annotation = []
    for i in range(len(ensg)):
        res = requests.get("http://mar2016.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r="%(ensg[i]))
        if "Synonyms" in res.text:
            soup = BeautifulSoup(res.text, 'html.parser')
            description = "<B>Description: </B>" + str(soup.find(class_='rhs').contents[0].contents[0]) + str(soup.find(class_='rhs').contents[0].contents[1]) + "]"
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[2].contents[0].contents[0].contents[0]) + " " + str(soup.find_all(class_='rhs')[2].contents[0].contents[1]) + " " + str(soup.find_all(class_='rhs')[2].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            synonym = "<B>Synonyms: </B>" + str(soup.find_all(class_='rhs')[1].contents[0].contents[0])
            link = "<a href = 'http://mar2016.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>"%(ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + description + "<br>" +  location + "<br>" + synonym + "<br>" + link + "<br></div>")
        elif "Description" not in res.text:
            soup = BeautifulSoup(res.text, 'html.parser')
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[0].contents[0].contents[0].contents[0]) + " " + str(
                soup.find_all(class_='rhs')[0].contents[0].contents[1]) + " " + str(
                soup.find_all(class_='rhs')[0].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            link = "<a href = 'http://mar2016.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>" % (
            ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + location  + "<br>" + link + "<br></div>")
        else:
            soup = BeautifulSoup(res.text, 'html.parser')
            description = "<B>Description: </B>" + str(soup.find(class_='rhs').contents[0].contents[0]) + str(
                soup.find(class_='rhs').contents[0].contents[1]) + "]"
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[1].contents[0].contents[0].contents[0]) + " " + str(
                soup.find_all(class_='rhs')[1].contents[0].contents[1]) + " " + str(
                soup.find_all(class_='rhs')[1].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            link = "<a href = 'http://mar2016.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>" % (
                ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + description + "<br>" + location + "<br>" + link + "<br></div>")
    final = []
    for i in annotation:
        final.append(Markup(i))
    return final

def mouse_annotations(ensg,genename):
    annotation = []
    for i in range(len(ensg)):
        res = requests.get("http://mar2016.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s;r="%(ensg[i]))
        if "Synonyms" in res.text:
            soup = BeautifulSoup(res.text, 'html.parser')
            description = "<B>Description: </B>" + str(soup.find(class_='rhs').contents[0].contents[0]) + str(soup.find(class_='rhs').contents[0].contents[1]) + "]"
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[2].contents[0].contents[0].contents[0]) + " " + str(soup.find_all(class_='rhs')[2].contents[0].contents[1]) + " " + str(soup.find_all(class_='rhs')[2].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            synonym = "<B>Synonyms: </B>" + str(soup.find_all(class_='rhs')[1].contents[0].contents[0])
            link = "<a href = 'http://mar2016.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>"%(ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + description + "<br>" +  location + "<br>" + synonym + "<br>" + link + "<br></div>")
        elif "Description" not in res.text:
            soup = BeautifulSoup(res.text, 'html.parser')
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[0].contents[0].contents[0].contents[0]) + " " + str(
                soup.find_all(class_='rhs')[0].contents[0].contents[1]) + " " + str(
                soup.find_all(class_='rhs')[0].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            link = "<a href = 'http://mar2016.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>" % (
            ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + location  + "<br>" + link+ "<br></div>")
        else:
            soup = BeautifulSoup(res.text, 'html.parser')
            description = "<B>Description: </B>" + str(soup.find(class_='rhs').contents[0].contents[0]) + str(
                soup.find(class_='rhs').contents[0].contents[1]) + "]"
            location = "<B>Location: </B>" + str(soup.find_all(class_='rhs')[1].contents[0].contents[0].contents[0]) + " " + str(
                soup.find_all(class_='rhs')[1].contents[0].contents[1]) + " " + str(
                soup.find_all(class_='rhs')[1].contents[1])
            location = location.replace("<p>", "<br>")
            location = location.replace("</p>", "")
            link = "<a href = 'http://mar2016.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r='>Go to ensembl database</a>" % (
                ensg[i])
            annotation.append("<div id='annotation_" + str(i) + "'>Annotations: <br>Gene Name: " + genename[i] + "<br>" + description + "<br>" + location + "<br>" + link + "<br></div>")
    final = []
    for i in annotation:
        final.append(Markup(i))
    return final

def transcript_expression(start, stop, chr, species):
    if species=="human":
        index = "human_transcript_exp"
        keys = ["Thyroid", "Testis", "Brain - Anterior cingulate cortex (BA24)", "Skin - Not Sun Exposed (Suprapubic)",
                "Esophagus - Mucosa", "Heart - Atrial Appendage", "Brain - Caudate (basal ganglia)",
                "Esophagus - Muscularis",
                "Brain - Putamen (basal ganglia)", "Small Intestine - Terminal Ileum", "Breast - Mammary Tissue",
                "Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Brain - Cerebellum", "Bladder",
                "Brain - Cerebellar Hemisphere", "Brain - Spinal cord (cervical c_1)", "Artery - Coronary", "Liver",
                "Esophagus - Gastroesophageal Junction", "Brain - Hypothalamus", "Colon - Transverse",
                "Brain - Amygdala",
                "Pancreas", "Adipose - Subcutaneous", "Cells - Leukemia cell line (CML)", "Spleen",
                "Brain - Hippocampus", "Whole Blood", "Brain - Cortex", "Artery - Tibial", "Uterus", "Stomach", "Ovary",
                "Artery - Aorta", "Heart - Left Ventricle", "Kidney - Cortex",
                "Brain - Nucleus accumbens (basalganglia)",
                "Prostate", "Brain - Frontal Cortex (BA9)", "Vagina", "Adipose - Visceral (Omentum)", "Adrenal Gland",
                "Lung",
                "Cells - Transformed fibroblasts", "Muscle - Skeletal", "Colon - Sigmoid", "Nerve - Tibial",
                "Brain - Substantia nigra", "Cells - EBV-transformed lymphocytes"]
        return_keys = ['Chr', 'Start', 'Stop', 'Transcript_stable_ID', 'Strand', 'Thyroid', 'Testis', 'Brain - Anterior cingulate cortex (BA24)', 'Skin - Not Sun Exposed (Suprapubic)', 'Esophagus - Mucosa', 'Heart - Atrial Appendage', 'Brain - Caudate (basal ganglia)', 'Esophagus - Muscularis', 'Brain - Putamen (basal ganglia)', 'Small Intestine - Terminal Ileum', 'Breast - Mammary Tissue', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Fallopian Tube', 'Brain - Cerebellum', 'Bladder', 'Brain - Cerebellar Hemisphere', 'Brain - Spinal cord (cervical c_1)', 'Artery - Coronary', 'Liver', 'Esophagus - Gastroesophageal Junction', 'Brain - Hypothalamus', 'Colon - Transverse', 'Brain - Amygdala', 'Strand', 'Pancreas', 'Adipose - Subcutaneous', 'Cells - Leukemia cell line (CML)', 'Spleen', 'Brain - Hippocampus', 'Whole Blood', 'Brain - Cortex', 'Artery - Tibial', 'Uterus', 'Stomach', 'Ovary', 'Artery - Aorta', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Brain - Nucleus accumbens (basalganglia)', 'Prostate', 'Brain - Frontal Cortex (BA9)', 'Vagina', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Lung', 'Cells - Transformed fibroblasts', 'Muscle - Skeletal', 'Colon - Sigmoid', 'Nerve - Tibial', 'Brain - Substantia nigra', 'Cells - EBV-transformed lymphocytes']
    else:
        index = "mouse_transcript_exp"
        keys = ["embryo", "heart", "bone marrow macrophage", "fat pad", "neural tube", "embryonic fibroblast", "brain", "hindbrain", "limb", "stomach", "erythroblast", "midbrain", "kidney", "B cell", "MEL cell line", "testis", "vesicular gland", "G1E", "subcutaneous adipose tissue", "adrenal gland", "gonadal fat pad", "telencephalon", "brown adipose tissue", "placenta", "intestine", "forestomach", "CH12.LX", "ES-Bruce4", "activated regulatory T-cells", "cortical plate", "regulatory T cell", "skeletal muscle tissue", "urinary bladder", "cerebellum", "small intestine", "416B", "NIH3T3", "pancreas", "A20", "Patski", "G1E-ER4", "embryonic facial prominence", "bone marrow", "spleen", "thymus", "splenic B cell", "inflammation-experienced regulatory T-cells", "forebrain", "uterus", "lung", "ovary", "muscle", "olfactory bulb", "liver"]
        return_keys = ['Chr','Start','Stop','Transcript_ID', 'Strand', 'embryo', 'heart', 'neural tube', 'bone marrow macrophage', 'CH12.LX', 'fat pad', 'embryonic fibroblast', 'brain', 'hindbrain', 'limb', 'stomach', 'erythroblast', 'kidney', 'B cell', 'MEL cell line', 'testis', 'vesicular gland', 'ES-Bruce4', 'G1E', 'subcutaneous adipose tissue', 'adrenal gland', 'gonadal fat pad', 'forestomach', 'brown adipose tissue', 'placenta', 'uterus', 'activated regulatory T-cells', 'intestine', 'cortical plate', 'regulatory T cell', 'skeletal muscle tissue', 'urinary bladder', 'embryonic facial prominence', 'small intestine', '416B', 'NIH3T3', 'midbrain', 'pancreas', 'cerebellum', 'Patski', 'G1E-ER4', 'bone marrow', 'spleen', 'thymus', 'A20', 'splenic B cell', 'telencephalon', 'forebrain', 'inflammation-experienced regulatory T-cells', 'lung', 'ovary', 'muscle', 'olfactory bulb', 'liver']

    request = []
    req_head = {'index': index, 'type': index}
    req_body = {
                "size": 9000,
                            "query": {
                                "bool": {
                                    "must": {
                                        "range": {"Start": {"gte": start, "lte": stop}}
                                    },
                                    "filter": {
                                        "term": {"Chr": chr}
                                    }
                                }
                            }
                        }
    request.extend([req_head, req_body])
    es = Elasticsearch('https://SECRET/elasticsearch', verify_certs=False, timeout=30,
                       max_retries=10, retry_on_timeout=True)
    resp = es.msearch(body=request)
    response = []
    expression_list = []
    try:
        for i in range(len(resp['responses'][0]['hits']['hits'])):
            response.append(resp['responses'][0]['hits']['hits'][i]["_source"])
            expression_list.append([])
            temp = resp['responses'][0]['hits']['hits'][i]["_source"]
            for j in keys:
                try:
                    expression_list[i].append(float(temp[j]))
                except:
                    expression_list[i].append(None)
    except:
        response = []
        expression_list = []
    normalized_list = []
    for i in range(len(expression_list)):
        normalized_list.append([])
        for j in expression_list[i]:
            if None not in expression_list[i] and max(expression_list[i])!=0:
                normalized_list[i].append(j/max(expression_list[i]))
            elif None in expression_list[i]:
                    maximum = 0
                    for k in expression_list[i]:
                        if k != None and k > maximum:
                            maximum = k
                    if maximum!=0 and j!=None:
                        normalized_list[i].append(j/maximum)
                    else:
                        normalized_list[i].append(j)
            else:
                normalized_list[i].append(j)
    return response, json.dumps(expression_list), json.dumps(normalized_list), return_keys

def get_genes(start,stop,chr):
    res = requests.post("https://SECRET/elasticsearch/annotations/_search?pretty=true&size=500",
                        verify=False, json=
                        {
                            "query": {
                                "bool": {
                                    "must": [{
                                        "range": {"START": {"gte": start, "lte": stop}},
                                        "range": {"STOP": {"lte": stop, "gte": start}}
                                    }],
                                    "filter": {
                                        "term": {"CHROM": chr}
                                    }
                                }
                            }
                        })
    data = res.json()["hits"]["hits"]
    genes = []
    ensembl_id = []
    for i in range(len(data)):
        if data[i]["_source"]["GENE_NAME"] not in genes:
            genes.append(data[i]["_source"]["GENE_NAME"])
            ensembl_id.append(data[i]["_source"]["GENE_ID"])
    return genes,ensembl_id

def get_genes_mouse(start,stop,chr):
    res = requests.post("https://SECRET/elasticsearch/mouse_annotations/_search?pretty=true&size=500",
                        verify=False, json=
                        {
                            "query": {
                                "bool": {
                                    "must": [{
                                        "range": {"Start": {"gte": start, "lte": stop}},
                                        "range": {"Stop": {"lte": stop, "gte": start}}
                                    }],
                                    "filter": {
                                        "term": {"Chr": chr}
                                    }
                                }
                            }
                        })
    data = res.json()["hits"]["hits"]
    genes = []
    ensembl_id = []
    for i in range(len(data)):
        if data[i]["_source"]["Gene"] not in genes:
            genes.append(data[i]["_source"]["Gene"])
            ensembl_id.append(data[i]["_source"]["ENSGID"])
    return genes,ensembl_id

def add_snp(mod, species):
    if mod!=[]:
        if species == "human":
            index = "snp"
            position = "POS"
        else:
            index = "mouse_snp"
            position = "START"
        chr = []
        start = []
        request = []
        for i in mod:
            if i!="!":
                start.append(i['Start'])
                chr.append(i['Chr'])
        for i in range(len(start)):
            req_head = {'index': index, 'type': index}
            req_body = {"size": 1,
                        "query": {
                            "bool": {
                                "must": [
                                    {
                                        "range": {
                                            position: {
                                                "gt": int(start[i])-1, "lt": int(start[i])+1
                                            }
                                        }
                                    },
                                    {
                                        "term": {
                                            "CHR": chr[i]
                                        }
                                    }]
                            }
                        }
                        }
            request.extend([req_head, req_body])
        es = Elasticsearch('https://SECRET/elasticsearch', verify_certs=False, timeout=50,
                           max_retries=10, retry_on_timeout=True)
        resp = es.msearch(body=request)
        response = {}
        for i in range(len(resp['responses'])):
            if resp['responses'][i]['hits']['hits'] != []:
                mod[i]["SNP"] = Markup("<a href='https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s' target='_blank'>%s</a>"%(str(resp['responses'][i]['hits']['hits'][0]['_source']['ID']).lstrip("rs"),resp['responses'][i]['hits']['hits'][0]['_source']['ID']))
            else:
                mod[i]["SNP"] = "-"
        return mod
    else:
        return mod

def mod_function(start, stop, chr, mod_indices, gene):
    request = []
    if gene!=0:
        for i in mod_indices:
            req_head = {'index': i, 'type': i}
            req_body = {
                        "size": 9000,
                          "query": {
                            "term" : { "Gene" : gene }
                          }
                        }
            request.extend([req_head, req_body])
    else:
        for i in mod_indices:
            req_head = {'index': i, 'type': i}
            req_body = {"size": 9000,
                                "query": {
                                    "bool": {
                                        "must": [{
                                            "range": {"Start": {"gte": start, "lte": stop}},
                                            "range": {"Stop": {"lte": stop, "gte": start}}
                                        }],
                                        "filter": {
                                            "term": {"Chr": chr}
                                        }
                                    }
                                }
                            }

            request.extend([req_head, req_body])
    es = Elasticsearch('https://SECRET/elasticsearch', verify_certs=False, timeout=50,
                       max_retries=10, retry_on_timeout=True)
    resp = es.msearch(body=request)
    response = {}
    for i in range(len(resp["responses"])):
        response[mod_indices[i]] = []
        for j in range(resp["responses"][i]["hits"]["total"]):
            try:
                resp["responses"][i]["hits"]["hits"][j]["_source"]["Gene"] = Markup("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['Gene'],resp['responses'][i]['hits']['hits'][j]['_source']['Gene']))
                resp["responses"][i]["hits"]["hits"][j]["_source"]["ENSG_ID"] = Markup("<a href='https://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r=' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['ENSG_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['ENSG_ID']))
                resp["responses"][i]["hits"]["hits"][j]["_source"]["ENST_ID"] = Markup("<a href='https://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s;r=' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['ENST_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['ENST_ID']))
                if "Hela" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"] = "HeLa"
                if "2360428323118480" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = "23604283"
                if "26863190" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = "26863196"
                if "/" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"] = resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"].replace("/", ", ")
                if "|" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]:
                    ids = resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]
                    ids = ids.split("|")
                    ids_final = []
                    for id in ids:
                        ids_final.append("<a href='https://www.ncbi.nlm.nih.gov/pubmed/%s' target='_blank'>%s</a>" % (id, id))
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = Markup(", ".join(ids_final))
                else:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = Markup("<a href='https://www.ncbi.nlm.nih.gov/pubmed/%s' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['Pubmed_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['Pubmed_ID']))
                response[mod_indices[i]].append(resp["responses"][i]["hits"]["hits"][j]["_source"])
            except IndexError:
                response[mod_indices[i]].append("!")
    pool_query = ThreadPool(processes=13)
    response['a_to_i_humans'] = pool_query.apply_async(add_snp, (response['a_to_i_humans'], 'human'))
    response['m1a_humans'] = pool_query.apply_async(add_snp, (response['m1a_humans'], 'human'))
    response['m5c_humans'] = pool_query.apply_async(add_snp, (response['m5c_humans'], 'human'))
    response['m6a_humans'] = pool_query.apply_async(add_snp, (response['m6a_humans'], 'human'))
    response['nm_humans'] = pool_query.apply_async(add_snp, (response['nm_humans'], 'human'))
    response['pseudou_humans'] = pool_query.apply_async(add_snp, (response['pseudou_humans'], 'human'))
    response['c_to_u_humans'] = pool_query.apply_async(add_snp, (response['c_to_u_humans'], 'human'))
    response['dihydrouridine_humans'] = pool_query.apply_async(add_snp, (response['dihydrouridine_humans'], 'human'))
    response['m1g_humans'] = pool_query.apply_async(add_snp, (response['m1g_humans'], 'human'))
    response['m2g_humans'] = pool_query.apply_async(add_snp, (response['m2g_humans'], 'human'))
    response['m7g_humans'] = pool_query.apply_async(add_snp, (response['m7g_humans'], 'human'))
    response['other_humans'] = pool_query.apply_async(add_snp, (response['other_humans'], 'human'))
    response['t6a_humans'] = pool_query.apply_async(add_snp, (response['t6a_humans'], 'human'))
    response['a_to_i_humans'] = response['a_to_i_humans'].get()
    response['m1a_humans'] = response['m1a_humans'].get()
    response['m5c_humans'] = response['m5c_humans'].get()
    response['m6a_humans'] = response['m6a_humans'].get()
    response['nm_humans'] = response['nm_humans'].get()
    response['pseudou_humans'] = response['pseudou_humans'].get()
    response['c_to_u_humans'] = response['c_to_u_humans'].get()
    response['dihydrouridine_humans'] = response['dihydrouridine_humans'].get()
    response['m1g_humans'] = response['m1g_humans'].get()
    response['m2g_humans'] = response['m2g_humans'].get()
    response['m7g_humans'] = response['m7g_humans'].get()
    response['other_humans'] = response['other_humans'].get()
    response['t6a_humans'] = response['t6a_humans'].get()
    pool_query.terminate()
    pool_query.close()
    return response['a_to_i_humans'],response['m1a_humans'],response['m5c_humans'],response['m6a_humans'],response['nm_humans'],response['pseudou_humans'], response['c_to_u_humans'], response['dihydrouridine_humans'], response['m1g_humans'],response['m2g_humans'], response['m7g_humans'], response['other_humans'], response['t6a_humans']

def mod_function_mouse(start, stop, chr, mod_indices, gene):
    request = []
    if gene!=0:
        for i in mod_indices:
            req_head = {'index': i, 'type': i}
            req_body = {
                        "size": 9000,
                          "query": {
                            "term" : { "Gene" : gene }
                          }
                        }
            request.extend([req_head, req_body])
    else:
        for i in mod_indices:
            req_head = {'index': i, 'type': i}
            req_body = {"size": 9000,
                                "query": {
                                    "bool": {
                                        "must": [{
                                            "range": {"Start": {"gte": start, "lte": stop}},
                                            "range": {"Stop": {"lte": stop, "gte": start}}
                                        }],
                                        "filter": {
                                            "term": {"Chr": chr}
                                        }
                                    }
                                }
                            }

            request.extend([req_head, req_body])
    es = Elasticsearch('https://SECRET/elasticsearch', verify_certs=False, timeout=30,
                       max_retries=10, retry_on_timeout=True)
    resp = es.msearch(body=request)
    response = {}
    for i in range(len(resp["responses"])):
        response[mod_indices[i]] = []
        for j in range(resp["responses"][i]["hits"]["total"]):
            try:
                resp["responses"][i]["hits"]["hits"][j]["_source"]["ENSG_ID"] = Markup("<a href='https://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s;r=' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['ENSG_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['ENSG_ID']))
                resp["responses"][i]["hits"]["hits"][j]["_source"]["ENST_ID"] = Markup("<a href='https://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s;r=' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['ENST_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['ENST_ID']))
                if "280771691" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = "28077169"
                if "/" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"]:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"] = resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"].replace("/", ", ")
                resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"] = resp["responses"][i]["hits"]["hits"][j]["_source"]["Tissue"].capitalize()

                if "|" in resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]:
                    ids = resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"]
                    ids = ids.split("|")
                    ids_final = []
                    for id in ids:
                        ids_final.append("<a href='https://www.ncbi.nlm.nih.gov/pubmed/%s' target='_blank'>%s</a>" % (id, id))
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = Markup(", ".join(ids_final))
                else:
                    resp["responses"][i]["hits"]["hits"][j]["_source"]["Pubmed_ID"] = Markup("<a href='https://www.ncbi.nlm.nih.gov/pubmed/%s' target='_blank'>%s</a>"%(resp['responses'][i]['hits']['hits'][j]['_source']['Pubmed_ID'],resp['responses'][i]['hits']['hits'][j]['_source']['Pubmed_ID']))
                response[mod_indices[i]].append(resp["responses"][i]["hits"]["hits"][j]["_source"])
            except IndexError:
                response[mod_indices[i]].append("!")
    pool_query = ThreadPool(processes=6)
    response['a_to_i_mouse'] = pool_query.apply_async(add_snp, (response['a_to_i_mouse'], 'mouse'))
    response['m1a_mouse'] = pool_query.apply_async(add_snp, (response['m1a_mouse'], 'mouse'))
    response['m5c_mouse'] = pool_query.apply_async(add_snp, (response['m5c_mouse'], 'mouse'))
    response['m6a_mouse'] = pool_query.apply_async(add_snp, (response['m6a_mouse'], 'mouse'))
    response['nm_mouse'] = pool_query.apply_async(add_snp, (response['nm_mouse'], 'mouse'))
    response['pseudou_mouse'] = pool_query.apply_async(add_snp, (response['pseudou_mouse'], 'mouse'))
    #response['c_to_u_mouse'] = pool_query.apply_async(add_snp, (response['c_to_u_mouse'], 'mouse'))
    response['dihydrouridine_mouse'] = pool_query.apply_async(add_snp, (response['dihydrouridine_mouse'], 'mouse'))
    response['m1g_mouse'] = pool_query.apply_async(add_snp, (response['m1g_mouse'], 'mouse'))
    response['m2g_mouse'] = pool_query.apply_async(add_snp, (response['m2g_mouse'], 'mouse'))
    response['m7g_mouse'] = pool_query.apply_async(add_snp, (response['m7g_mouse'], 'mouse'))
    response['other_mouse'] = pool_query.apply_async(add_snp, (response['other_mouse'], 'mouse'))
    response['t6a_mouse'] = pool_query.apply_async(add_snp, (response['t6a_mouse'], 'mouse'))
    response['a_to_i_mouse'] = response['a_to_i_mouse'].get()
    response['m1a_mouse'] = response['m1a_mouse'].get()
    response['m5c_mouse'] = response['m5c_mouse'].get()
    response['m6a_mouse'] = response['m6a_mouse'].get()
    response['nm_mouse'] = response['nm_mouse'].get()
    response['pseudou_mouse'] = response['pseudou_mouse'].get()
    #response['c_to_u_mouse'] = response['c_to_u_mouse'].get()
    response['dihydrouridine_mouse'] = response['dihydrouridine_mouse'].get()
    response['m1g_mouse'] = response['m1g_mouse'].get()
    response['m2g_mouse'] = response['m2g_mouse'].get()
    response['m7g_mouse'] = response['m7g_mouse'].get()
    response['other_mouse'] = response['other_mouse'].get()
    response['t6a_mouse'] = response['t6a_mouse'].get()
    pool_query.terminate()
    pool_query.close()
    return response['a_to_i_mouse'],response['m1a_mouse'],response['m5c_mouse'],response['m6a_mouse'],response['nm_mouse'],response['pseudou_mouse'], response['dihydrouridine_mouse'], response['m1g_mouse'],response['m2g_mouse'], response['m7g_mouse'], response['other_mouse'], response['t6a_mouse']

app = Flask(__name__)

#___________________________________________________API SECTION STARTS HERE____________________________________________________________________________________#
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = True

api_base_path = '/home/sliceit/epitomy/static/api_data'
#api_base_path = "C:\\Users\\MeanMachine\\Desktop\\rna_mods"


@app.route('/api/<string:species>/<string:modification>/position/<int:position>', methods=['GET'])
def pos(species,modification,position):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    start = [task for task in tasks if task['Start'] == position]
    if len(start) == 0:
        abort(404)
    return jsonify({'Position': start})

@app.route('/api/<string:species>/<string:modification>/chromosome/<string:chr>', methods=['GET'])
def chromosome(species, modification, chr):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    chrom = [task for task in tasks if str(task['Chr']).lower() == str(chr).lower()]
    if len(chrom) == 0:
        abort(404)
    return jsonify({'Chromosome': chrom})

@app.route('/api/<string:species>/<string:modification>/tissue/<string:tissue>', methods=['GET'])
def tissues(species, modification, tissue):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    Tissue = [task for task in tasks if str(tissue).lower() in str(task['Tissue']).lower()]
    if len(Tissue) == 0:
        abort(404)
    return jsonify({'Tissue': Tissue})

@app.route('/api/<string:species>/<string:modification>/gene/<string:gene>', methods=['GET'])
def genes(species, modification, gene):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    Gene = [task for task in tasks if str(gene).lower() in str(task['Gene']).lower()]
    if len(Gene) == 0:
        abort(404)
    return jsonify({'Gene': Gene})

@app.route('/api/<string:species>/<string:modification>/pubmed/<string:pubmedid>', methods=['GET'])
def pubmed(species, modification, pubmedid):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    Pubmed = [task for task in tasks if str(pubmedid).lower() in str(task['PubmedID']).lower()]
    if len(Pubmed) == 0:
        abort(404)
    return jsonify({'Pubmed': Pubmed})

@app.route('/api/<string:species>/<string:modification>/multi/<string:chr>_<string:gene>_<string:tissue>_<string:pubmedid>', methods=['GET'])
def multi(species, modification, gene, chr, tissue, pubmedid):
    with open('{}/{}/{}.json'.format(api_base_path, species, modification)) as json_data:
        tasks = json.load(json_data)
    if gene == "*":
        gene = ""
    if chr == "*":
        chr = ""
    if tissue == "*":
        tissue = ""
    if pubmedid == "*":
        pubmedid = ""
    if chr != "":
        out = [task for task in tasks if str(gene).lower() in str(task['Gene']).lower() and str(chr).lower() == task['Chr'] and str(tissue).lower() in str(task['Tissue']).lower() and str(pubmedid).lower() in str(task['PubmedID']).lower()]
    else:
        out = [task for task in tasks if str(gene).lower() in str(task['Gene']).lower() and str(tissue).lower() in str(task['Tissue']).lower() and str(pubmedid).lower() in str(task['PubmedID']).lower()]
    if len(out) == 0:
        abort(404)
    return jsonify(out)

@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify([{'Chr': 'Not found', 'Start': 'Not found', 'Strand': 'Not found', 'Gene': 'Not found', 'Tissue': 'Not found', 'PubmedID': 'Not found'}]))

#___________________________________________________API SECTION ENDS HERE____________________________________________________________________________________#


@app.route('/contact', methods=['GET', 'POST'])
def contact():
    if request.method == "POST":
        fname = request.form["fname"]
        lname = request.form["lname"]
        femail = request.form["email"]
        reason = request.form["contact_reason"]
        msg = str(request.form["msg"]).encode("utf-8")
        website = request.form["website"]
        body = """
        Website: %s 
        First Name: %s 
        Last Name %s 
        Email: %s 
        Contact reason: %s 
        Message: %s"""%(website, fname, lname, femail, reason, msg)

        email = "jangalab.iupui@gmail.com"
        pwd = "Hesoyam1"

        smtp_server = smtplib.SMTP('smtp.gmail.com', 587)
        smtp_server.ehlo()
        smtp_server.starttls()
        smtp_server.login(email, pwd)

        to_emails = ["savemuri@iu.edu", "scjanga@iupui.edu"]

        for i in to_emails:
            smtp_server.sendmail(email, i, body)

        smtp_server.quit()
        alert = Markup('<div class="container"><div class="fire_alert alert alert-success alert-dismissable"><button aria-hidden="true" data-dismiss="alert" class="close" type="button">×</button><strong>Message has been successfully sent!</strong></div></div>')
        return render_template('contact.html', alert=alert)

    return render_template('contact.html', alert=0)

@app.route('/user_guide', methods=['GET', 'POST'])
def user_guide():
    return render_template('user_guide.html')

@app.route('/documentation', methods=['GET', 'POST'])
def documentation():
    return render_template('documentation.html')

@app.route('/', methods=['POST', 'GET'])
def process_data():
    human_indices = ["a_to_i_humans", "m1a_humans", "m5c_humans", "m6a_humans", "nm_humans",
                   "pseudou_humans", "c_to_u_humans", "dihydrouridine_humans", "m1g_humans", "m2g_humans", "m7g_humans", "other_humans", "t6a_humans"]
    mouse_indices = ["a_to_i_mouse", "m1a_mouse", "m5c_mouse", "m6a_mouse", "nm_mouse",
                     "pseudou_mouse", "dihydrouridine_mouse", "m1g_mouse", "m2g_mouse", "m7g_mouse", "other_mouse", "t6a_mouse"]
    file_names = ["A_to_I", "m1a", "m5c", "m6a", "nm", "pseudou"]
    all_mods = ['a_to_i', 'm1a', 'm5c', 'm6a', 'nm', 'pseudou', 'dihydrouridine', 'm1g', 'm2g', 'm7g', 'other', 't6a']

    modification_headers = ["Chromosome", "Modification Co-ordinate", "SNP", "Tissue source", "PubMed ID", "Strand", "ENSEMBL Gene ID", "ENSEMBL Transcript ID", "Gene symbol", "RNA Type", "Base before modification"]
    if request.method == "POST":
        if request.form.get("check") == "human_gene":
            try:
                species = "human"
                gene = str(request.form["gene_search_field"])
                gene, ensembl_id = gene.split(" - ")[0], gene.split(" - ")[1]
                gene = gene.lower()
                annotations = human_annotations([ensembl_id], [gene.upper()])
                browser_coordinates = genes_coordinates[gene.upper()]
                chr = browser_coordinates.split(":")[0]
                start = browser_coordinates.split(":")[1].split("-")[0]
                stop = browser_coordinates.split(":")[1].split("-")[1]
                location_coordinates = "%s:%s-%s"%(chr,start,stop)
                a_to_i, m1a, m5c, m6a, nm, pseudou, c_to_u, dihydrouridine, m1g, m2g, m7g, other, t6a  = mod_function(start,stop,chr,human_indices,gene)
                expression_final, expression_list, normalized_list, expression_keys = transcript_expression(start,stop,chr,species)
                expression_files = ["Thyroid","Testis","Brain-AnteriorcingulatecortexBA24","Skin-NotSunExposedSuprapubic",
                                    "Esophagus-Mucosa","Heart-AtrialAppendage","Brain-Caudatebasalganglia",
                                    "Esophagus-Muscularis",
                                    "Brain-Putamenbasalganglia","SmallIntestine-TerminalIleum","Breast-MammaryTissue",
                                    "Cervix-Ectocervix","Cervix-Endocervix","FallopianTube","Brain-Cerebellum","Bladder",
                                    "Brain-CerebellarHemisphere","Brain-Spinalcordcervicalc_1","Artery-Coronary","Liver",
                                    "Esophagus-GastroesophagealJunction","Brain-Hypothalamus","Colon-Transverse",
                                    "Brain-Amygdala",
                                    "Pancreas","Adipose-Subcutaneous","Cells-LeukemiacelllineCML","Spleen",
                                    "Brain-Hippocampus","WholeBlood","Brain-Cortex","Artery-Tibial","Uterus","Stomach","Ovary",
                                    "Artery-Aorta","Heart-LeftVentricle","Kidney-Cortex",
                                    "Brain-Nucleusaccumbensbasalganglia",
                                    "Prostate","Brain-FrontalCortexBA9","Vagina","Adipose-VisceralOmentum","AdrenalGland",
                                    "Lung",
                                    "Cells-Transformedfibroblasts","Muscle-Skeletal","Colon-Sigmoid","Nerve-Tibial",
                                    "Brain-Substantianigra","Cells-EBV-transformedlymphocytes"]
                return render_template("result.html", annotations=annotations, a_to_i=a_to_i, m1a=m1a, m5c=m5c, m6a=m6a, nm=nm, pseudou=pseudou, c_to_u=c_to_u, dihydrouridine=dihydrouridine, m1g=m1g, m2g=m2g, m7g=m7g, other=other, t6a=t6a, modification_headers=modification_headers, location_coordinates=location_coordinates, species=species, expression_final=expression_final, expression_list=expression_list, normalized_list=normalized_list, file_names=file_names, expression_keys=expression_keys, expression_files=expression_files, all_mods=all_mods)
            except:
                return render_template("error.html")
        elif request.form.get("check") == "human_coordinates":
            try:
                species = "human"
                coordinates = request.form["gene_search_field"]
                chr = request.form.get("chrom")
                start = coordinates.split("|")[0]
                stop = coordinates.split("|")[1]
                location_coordinates = "%s:%s-%s"%(chr,start,stop)
                if int(stop)-int(start) >100000:
                    stop = int(start)+100000
                genes, ensembl_id = get_genes(start,stop,chr)
                annotations = human_annotations(ensembl_id, genes)
                a_to_i, m1a, m5c, m6a, nm, pseudou, c_to_u, dihydrouridine, m1g, m2g, m7g, other, t6a = mod_function(start,stop,chr,human_indices, 0)
                expression_final, expression_list, normalized_list, expression_keys = transcript_expression(start,stop,chr,species)
                expression_files = ["Thyroid","Testis","Brain-AnteriorcingulatecortexBA24","Skin-NotSunExposedSuprapubic",
                                    "Esophagus-Mucosa","Heart-AtrialAppendage","Brain-Caudatebasalganglia",
                                    "Esophagus-Muscularis",
                                    "Brain-Putamenbasalganglia","SmallIntestine-TerminalIleum","Breast-MammaryTissue",
                                    "Cervix-Ectocervix","Cervix-Endocervix","FallopianTube","Brain-Cerebellum","Bladder",
                                    "Brain-CerebellarHemisphere","Brain-Spinalcordcervicalc_1","Artery-Coronary","Liver",
                                    "Esophagus-GastroesophagealJunction","Brain-Hypothalamus","Colon-Transverse",
                                    "Brain-Amygdala",
                                    "Pancreas","Adipose-Subcutaneous","Cells-LeukemiacelllineCML","Spleen",
                                    "Brain-Hippocampus","WholeBlood","Brain-Cortex","Artery-Tibial","Uterus","Stomach","Ovary",
                                    "Artery-Aorta","Heart-LeftVentricle","Kidney-Cortex",
                                    "Brain-Nucleusaccumbensbasalganglia",
                                    "Prostate","Brain-FrontalCortexBA9","Vagina","Adipose-VisceralOmentum","AdrenalGland",
                                    "Lung",
                                    "Cells-Transformedfibroblasts","Muscle-Skeletal","Colon-Sigmoid","Nerve-Tibial",
                                    "Brain-Substantianigra","Cells-EBV-transformedlymphocytes"]
                return render_template("result.html",  a_to_i=a_to_i, m1a=m1a, m5c=m5c, m6a=m6a, nm=nm, pseudou=pseudou, c_to_u=c_to_u, dihydrouridine=dihydrouridine, m1g=m1g, m2g=m2g, m7g=m7g, other=other, t6a=t6a, modification_headers=modification_headers, annotations=annotations, location_coordinates=location_coordinates, species=species, expression_final=expression_final, expression_list=expression_list, normalized_list=normalized_list, file_names=file_names, expression_keys=expression_keys, expression_files=expression_files, all_mods=all_mods)
            except:
                return render_template("error.html")
        elif request.form.get("check") == "mouse_gene":
            try:
                species = "mouse"
                gene = str(request.form["gene_search_field"])
                gene, ensembl_id = gene.split(" - ")[0], gene.split(" - ")[1]
                gene = gene.lower()
                annotations = mouse_annotations([ensembl_id], [gene.lower()])
                browser_coordinates = genes_coordinates_mouse[gene.lower()]
                chr = browser_coordinates.split(":")[0]
                chr_out = chr
                if chr == 'chrx':
                    chr_out = "chrX"
                elif chr == 'chry':
                    chr_out = "chrY"
                start = browser_coordinates.split(":")[1].split("-")[0]
                stop = browser_coordinates.split(":")[1].split("-")[1]
                location_coordinates = "%s:%s-%s"%(chr_out,start,stop)
                a_to_i, m1a, m5c, m6a, nm, pseudou, dihydrouridine, m1g, m2g, m7g, other, t6a = mod_function_mouse(start,stop,chr,mouse_indices,gene)
                expression_final, expression_list, normalized_list, expression_keys = transcript_expression(start,stop,chr,species)
                expression_files = ['embryo', 'heart', 'bonemarrowmacrophage', 'fatpad', 'neuraltube', 'embryonicfibroblast', 'brain', 'hindbrain', 'limb', 'stomach', 'erythroblast', 'midbrain', 'kidney', 'Bcell', 'MELcellline', 'testis', 'vesiculargland', 'G1E', 'subcutaneousadiposetissue', 'adrenalgland', 'gonadalfatpad', 'telencephalon', 'brownadiposetissue', 'placenta', 'intestine', 'forestomach', 'CH12.LX', 'ES-Bruce4', 'activatedregulatoryT-cells', 'corticalplate', 'regulatoryTcell', 'skeletalmuscletissue', 'urinarybladder', 'cerebellum', 'smallintestine', '416B', 'NIH3T3', 'pancreas', 'A20', 'Patski', 'G1E-ER4', 'embryonicfacialprominence', 'bonemarrow', 'spleen', 'thymus', 'splenicBcell', 'inflammation-experiencedregulatoryT-cells', 'forebrain', 'uterus', 'lung', 'ovary', 'muscle', 'olfactorybulb', 'liver']
                return render_template("result.html", annotations=annotations, a_to_i=a_to_i, m1a=m1a, m5c=m5c, m6a=m6a, nm=nm, pseudou=pseudou, dihydrouridine=dihydrouridine, m1g=m1g, m2g=m2g, m7g=m7g, other=other, t6a=t6a, modification_headers=modification_headers, location_coordinates=location_coordinates, species=species, expression_final=expression_final, expression_list=expression_list, normalized_list=normalized_list, file_names=file_names, expression_keys=expression_keys, expression_files=expression_files, all_mods=all_mods)
            except:
                return render_template("error.html")
        elif request.form.get("check") == "mouse_coordinates":
            try:
                species = "mouse"
                coordinates = request.form["gene_search_field"]
                chr = request.form.get("chrom")
                chr_out = chr
                if chr == 'chrx':
                    chr_out = "chrX"
                elif chr == 'chry':
                    chr_out = "chrY"
                start = coordinates.split("|")[0]
                stop = coordinates.split("|")[1]
                location_coordinates = "%s:%s-%s"%(chr_out,start,stop)
                if int(stop)-int(start) >100000:
                    stop = int(start)+100000
                genes, ensembl_id = get_genes_mouse(start,stop,chr)
                annotations = mouse_annotations(ensembl_id, genes)
                a_to_i, m1a, m5c, m6a, nm, pseudou, dihydrouridine, m1g, m2g, m7g, other, t6a = mod_function_mouse(start,stop,chr,mouse_indices,0)
                expression_final, expression_list, normalized_list, expression_keys = transcript_expression(start,stop,chr,species)
                expression_files = ['embryo', 'heart', 'bonemarrowmacrophage', 'fatpad', 'neuraltube', 'embryonicfibroblast', 'brain', 'hindbrain', 'limb', 'stomach', 'erythroblast', 'midbrain', 'kidney', 'Bcell', 'MELcellline', 'testis', 'vesiculargland', 'G1E', 'subcutaneousadiposetissue', 'adrenalgland', 'gonadalfatpad', 'telencephalon', 'brownadiposetissue', 'placenta', 'intestine', 'forestomach', 'CH12.LX', 'ES-Bruce4', 'activatedregulatoryT-cells', 'corticalplate', 'regulatoryTcell', 'skeletalmuscletissue', 'urinarybladder', 'cerebellum', 'smallintestine', '416B', 'NIH3T3', 'pancreas', 'A20', 'Patski', 'G1E-ER4', 'embryonicfacialprominence', 'bonemarrow', 'spleen', 'thymus', 'splenicBcell', 'inflammation-experiencedregulatoryT-cells', 'forebrain', 'uterus', 'lung', 'ovary', 'muscle', 'olfactorybulb', 'liver']
                return render_template("result.html",  a_to_i=a_to_i, m1a=m1a, m5c=m5c, m6a=m6a, nm=nm, pseudou=pseudou, dihydrouridine=dihydrouridine, m1g=m1g, m2g=m2g, m7g=m7g, other=other, t6a=t6a, modification_headers=modification_headers, annotations=annotations, location_coordinates=location_coordinates, species=species, expression_final=expression_final, expression_list=expression_list, normalized_list=normalized_list, file_names=file_names, expression_keys=expression_keys, expression_files=expression_files, all_mods=all_mods)
            except:
                return render_template("error.html")
        elif request.form.get("check") == "mod_type" or request.form.get("check") == "mouse_mod_type" or request.form.get("check") == "pubmed_type" or request.form.get("check") == "mouse_pubmed_type":
                if request.form.get("check") == "mouse_mod_type" or request.form.get("check") == "mouse_pubmed_type":
                    species = "mouse"
                elif request.form.get("check") == "mod_type" or request.form.get("check") == "pubmed_type":
                    species = "human"
                modification = request.form.get("modification")
                chr = "*"
                gene = "*"
                tissue = "*"
                if "|" in modification:
                    modification, pubmedid = modification.split("|")
                else:
                    pubmedid = "*"
                # CHR - GENE - TISSUE - PUBMEDID
                url = "https://epitomy.soic.iupui.edu/api/{}/{}/multi/{}_{}_{}_{}".format(species, modification, chr, gene, tissue, pubmedid)
                #url = "http://127.0.0.1:5000/api/{}/{}/multi/{}_{}_{}_{}".format(species, modification, chr, gene, tissue, pubmedid)
                return render_template("mod_search.html", url=url, modification=modification)

    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)