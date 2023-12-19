import streamlit as st
import pandas as pd
import numpy as np
import pandas.io.sql as sqlio
from utils_streamlit import run_optics_table

from io import StringIO
import json

from database import *
from sila2.client import SilaClient
import time
import psycopg2 as psycopg
from pathlib import Path
import os
import sila2
from sila2.client import SilaClient
import time
import sys
import json
from database import *
from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    ForeignKey,
    DateTime,
    LargeBinary

)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSON
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.functions import  current_timestamp
from datetime import datetime
import json
import streamlit as st
import psycopg2
from pyunpack import Archive
from shutil import move
from PIL import Image


CHEMSPEED_IP = "127.0.0.1"
CHEMSPEED_CERT = None
CHEMSPEED_PORT = 65001
HPLC_IP = "127.0.0.1"
HPLC_CERT = None
HPLC_PORT = 65010
OPTICS_IP = "127.0.0.1"
OPTICS_CERT = None
OPTICS_PORT = 65070
ATLAS_IP = "127.0.0.1"
ATLAS_PORT = 65100


from dblogin import get_login
dbname, dbuser, dbpassword = get_login()

engine = create_engine(f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}")
database_connect = f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}"




tab1, tab2 = st.tabs(["Run Job", "analyze data"])


@st.cache_resource
def init_connection():
    return psycopg2.connect(f"dbname={dbname} user={dbuser} password={dbpassword}")

conn = init_connection()

@st.cache_data(ttl=1)
def get_jobs(device, date):
    if len(device) == 1:
        device, = device
        with conn.cursor() as cur:
            cur.execute(f"SELECT name, id, timestamp FROM job WHERE devices_id={device} AND timestamp BETWEEN '{date}T00:00:00.00' AND '{date}T23:59:59.999'")
            return cur.fetchall()   
    else:
        with conn.cursor() as cur:
            cur.execute(f"SELECT name, id, timestamp FROM job WHERE devices_id IN {device} AND timestamp BETWEEN '{date}T00:00:00.00' AND '{date}T23:59:59.999'")
            return cur.fetchall()


@st.cache_data(ttl=1)
def fetch_optics_results(job_id):

    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT code FROM opticsdata WHERE job_id={job_id}")
            data,  = cur.fetchall()[0]
        
        if not os.path.exists(f"pages/bigjobresults/optics/{job_id}"):
            os.mkdir(f"pages/bigjobresults/optics/{job_id}")
        
        if data is None:
            return False
        
        with open(f"pages/bigjobresults/optics/{job_id}/{job_id}.7z", "wb") as f:
            f.write(data)
        
        if os.path.exists(f"pages/bigjobresults/optics/{job_id}/extractedconfirm.txt"):
            return True
        
        else:
            os.system(f"7za e pages/bigjobresults/optics/{job_id}/{job_id}.7z  -opages/bigjobresults/optics/{job_id}")

            with open(f"pages/bigjobresults/optics/{job_id}/extractedconfirm.txt", "w") as f:
                f.write("file extracted")

        # root = f'pages/bigjobresults/optics/{job_id}/output_folder'
        # child = os.listdir(root)[0]
        # try:
        #     os.rmdir(os.path.join(root, child))
        # except:
        #     return True

        return True


    except:
        return False
    


with tab1:
    st.title("Optical Table Interface")
    with st.form("OpticsJob"):
        st.subheader("Job Submission to Optics Table")

        st.write("Please name the job in the job file")

        
        uploader = st.file_uploader("upload your job file", type="json")

        placeholder = st.empty()
        
        with placeholder.container():
            submitted = st.form_submit_button("Submit")
        
        if submitted:

            data = json.load(uploader)

            jobname = data["name"]

            jobtimestamp = datetime.now()
            engine = create_engine(database_connect)
            Session = sessionmaker(bind=engine)
            session = Session()


            with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbpassword}") as conn:
                with conn.cursor() as cur:
                    cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
                    VALUES (
                        '{data["name"]}',
                        '{3}',
                        '{json.dumps(data)}',
                        '{"test hplc job"}',
                        '{jobtimestamp}',
                        '{"uoft"}'    
                        );
                    """)
                    conn.commit() 

            JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

            st.write(f"your job id is {JOB_ID}")

            instance = run_optics_table(data, JOB_ID, OPTICS_IP, OPTICS_PORT, OPTICS_CERT)


            st.warning("Job submitted successfully- to view results: go to display data. to submit a new job, refresh the page or click on the button below")

            placeholder.empty()

            submitted = st.form_submit_button("re-fresh and submit new job")

            if submitted:
                st.experimental_rerun()

with tab2:
    st.title("Optical Table Interface")
    with st.form("findjob"):
        st.subheader("Find jobs")
        df = sqlio.read_sql_query("select id, name from device", conn)
        date = st.date_input("choose date of job")

        submitted = st.form_submit_button("Submit")
        if submitted:
            data = st.dataframe(pd.DataFrame(get_jobs(tuple([3]), date), columns=["name", "job id", "time"]))


    st.subheader("Enter job id")
   
    job_id = st.number_input("job_id- must be an integer!!", min_value=0, value=0, step=1)

    if job_id == 0:
        st.write(f"return data not found for job {job_id}")
    else:
        filefound = fetch_optics_results(job_id)

        if filefound:

            root = f"pages/bigjobresults/optics/{job_id}"

            st.write(f"return data found for job {job_id}")
            with open(f"{root}/{job_id}.7z", "rb") as f:
                data = f.read()
            st.download_button("download compressed data", data)

            files = os.listdir(root)

            for f in files:
                if f.endswith('.png'):

                    image = Image.open(os.path.join(root,f))
                    st.image(image)

        else:
            st.write(f"return data not found for job {job_id}")

            





    



        

        


