import streamlit as st
from utils_streamlit import hplc_and_chemspeed_inject
import json
from database import Job
import psycopg2 as psycopg
import json
from sqlalchemy import (
    create_engine,
)
import time
import pandas as pd
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import streamlit as st
import psycopg2
import os
from pathlib import Path
import pickle
from PIL import Image
from dblogin import get_login

dbname, dbuser, dbpassword = get_login()
engine = create_engine(f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}")
database_connect = f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}"


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


#@st.cache_data(ttl=1, persist=True)
def fetch_hplc_results(job_id):
    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT code FROM hplcdata WHERE job_id={job_id}")
            data,  = cur.fetchall()[0]
        
        if not os.path.exists(f"pages/bigjobresults/hplc/{job_id}"):
            os.mkdir(f"pages/bigjobresults/hplc/{job_id}")
        
        if data is None:
            return False
         
        with open(f"pages/bigjobresults/hplc/{job_id}/{job_id}.7z", "wb") as f:
            f.write(data)

        confirmfile = Path(f"pages/bigjobresults/hplc/{job_id}/confirm")

        if confirmfile.is_dir():
            return True

        else:
            os.system(f"7za e pages/bigjobresults/hplc/{job_id}/{job_id}.7z -aos -opages/bigjobresults/hplc/{job_id}")
            os.mkdir(f"pages/bigjobresults/hplc/{job_id}/confirm")   
        return True
    except:
       return False



with tab1:
    st.title("HPLC interface")
    with st.form("HPLCJob"):
        st.subheader("Job Submission to HPLC")

        st.write("To execute an injection from chemspeed to the HPLC, first create an hplc job below. Then, wait for the HPLC status to be ready. enter the position to inject from chemspeed")

        jobname = st.text_input("name your job", max_chars=100)
        uploader = st.file_uploader("upload your job file", type="json")
        injectposition = st.text_input("position on chemspeed to inject", max_chars=7)

        placeholder = st.empty()
        
        with placeholder.container():
            submitted = st.form_submit_button("Submit")
        
        if submitted:
            if jobname == None or uploader == None:
                st.text("jobname is invalid")
            
            else:

                st.write("please wait, submitting job..")

                data = json.load(uploader)

                jobtimestamp = datetime.now()
                engine = create_engine(database_connect)
                Session = sessionmaker(bind=engine)
                session = Session()


               
                hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(data, injectposition)
                
                st.write(f"job id is: {resultsdict}")

                st.warning("Job submitted successfully- to view results: go to display data. to submit a new job, refresh the page or click on the button below")

                

                while chemspeed_thread.is_alive():
                    time.sleep(1)

                placeholder.empty()

                submitted = st.form_submit_button("re-fresh and submit new job")

                if submitted:
                    st.experimental_rerun()


with tab2:

    st.title("HPLC Table Interface")

    st.subheader("Scroll through table to see jobs")
    with conn.cursor() as cur:
        cur.execute(f"SELECT name, id, timestamp, code FROM job WHERE devices_id={2}")
        df = pd.DataFrame(cur.fetchall(), columns=["name", "job id", "time", "code"])  

    st.dataframe(df)



    st.subheader("Enter job id and Visualize results")
   
    job_id = st.number_input("job_id- must be an integer!!", min_value=0, value=0, step=1)

    if job_id == 0:
        st.write(f"return data not found for job {job_id}")
    else:
        filefound = fetch_hplc_results(job_id)

        if filefound:

            root = f"pages/bigjobresults/hplc/{job_id}"
            #outputfolder = os.path.join(root, "output_folder")

            st.write(f"return data found for job {job_id}")
            with open(f"{root}/{job_id}.7z", "rb") as f:
                data = f.read()
            st.download_button("download compressed data", data)

            files = os.listdir(root)

            for f in files:
                if f.endswith('.png'):
                    image = Image.open(os.path.join(root,f))
                    st.image(image)
                if f.endswith('.csv'):
                    data = st.dataframe(pd.read_csv(os.path.join(root,f)))
                if f== "characterization_params.json":
                    with open(os.path.join(root,f), "r") as f:
                        jobfile = json.load(f)
                    st.write("Optics job file:")
                    st.json(jobfile, expanded = True)
                if f=="target_retention_times.pkl":
                    with open(os.path.join(root,f), "rb") as f:
                        jobfile = pickle.load(f)
                    st.write("retention times:")
                    st.json(jobfile, expanded = True)


        else:
            st.write(f"return data not found for job {job_id}")

            
    






