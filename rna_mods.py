import os
os.chdir("/home/sliceit/epitomy/rna_mods/hg38/")

import mysql.connector
from sqlalchemy import create_engine

user = 'sliceit_epitomy'
passw = 'Hesoyam1!'
host =  'localhost'
port = 3306
database = 'sliceit_epitomy'

mydb = create_engine('mysql+pymysql://' + user + ':' + passw + '@' + host + ':' + str(port) + '/' + database , echo=False)

a = os.listdir("/home/sliceit/epitomy/rna_mods/hg38/")

import pandas as pd

for i in a:
    df = pd.read_csv("%s"%i, sep="\t", )
    del df['Stop']
    df.to_sql(con=mydb, name=str(i).rstrip(".txt"), if_exists='replace', index=False)


os.chdir("/home/sliceit/epitomy/rna_mods/mm10/")

mydb = create_engine('mysql+pymysql://' + user + ':' + passw + '@' + host + ':' + str(port) + '/' + database , echo=False)

a = os.listdir("/home/sliceit/epitomy/rna_mods/mm10/")

for i in a:
    df = pd.read_csv("%s"%i, sep="\t", )
    del df['Stop']
    df.to_sql(con=mydb, name=str(i).rstrip(".txt"), if_exists='replace', index=False)
