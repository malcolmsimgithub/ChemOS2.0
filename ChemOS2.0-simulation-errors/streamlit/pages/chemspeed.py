import streamlit as st
import pandas as pd
import numpy as np

from streamlitcommands import run_chemspeed_synthesis

from io import StringIO
import json

from database import *
from sila2.client import SilaClient
import time
import psycopg2 as psycopg
import pathlib
from pathlib import Path
import os
import uuid
import sila2
from sila2.client import SilaClient
import time
import sys
import json
from database import *
from sqlalchemy import (
    create_engine,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSON
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import threading
import streamlit as st
import threading

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

st.title("Chemspeed")

        

with st.form("ChemspeedJob"):
    st.subheader("Batch Synthesis")
    st.write("Job Submission to chemspeed")

    jobname = st.text_input("name your job", max_chars=100)
    description = st.text_input("job description", max_chars= 1000)
    uploader = st.file_uploader("upload your job file", type="json")

    placeholder = st.empty()
    
    with placeholder.container():
        submitted = st.form_submit_button("Submit")
    
    if submitted:
        if jobname == None or uploader == None:
            st.text("jobname is invalid")
        
        else:

            data = json.load(uploader)

            jobtimestamp = datetime.now()
            engine = create_engine(database_connect)
            Session = sessionmaker(bind=engine)
            session = Session()

            with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbpassword}") as conn:
                with conn.cursor() as cur:
                    cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
                    VALUES (
                        '{jobname}',
                        '{1}',
                        '{json.dumps(data)}',
                        '{description}',
                        '{jobtimestamp}',
                        '{"uoft"}'    
                        );
                    """)
                conn.commit() 
            
            JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

            st.write(f"your job id is {JOB_ID}")

            instance = run_chemspeed_synthesis(data, jobname, JOB_ID, CHEMSPEED_IP, CHEMSPEED_PORT, CHEMSPEED_CERT)


            st.warning("Job submitted successfully- to view results: go to display data. to submit a new job, refresh the page or click on the button below")

            placeholder.empty()

            submitted = st.form_submit_button("re-fresh and submit new job")

            if submitted:
                st.experimental_rerun()

with st.form("Get Inventory and Rack Positions"):
    st.subheader("See the inventory and storage rack positions of the chemspeed")
   
    submitted = st.form_submit_button("Show")
    
    if submitted:
        try:
            st.write("click on arrows to expand")
            client = SilaClient(CHEMSPEED_IP, 65002, root_certs=open(CHEMSPEED_CERT, "rb").read())

            response = client.ChemSpeedOperator.GetRackpositions()

            dict = eval(client.ChemSpeedOperator.GetInventory().Termination)

            st.write("Inventory:")
            st.json(dict, expanded = False)

            st.write("RACKL:")
            st.json(eval(response.RACKL), expanded= False)

            st.write("RACKR:")
            st.json(eval(response.RACKR), expanded = False)

            st.write("SPE:")
            st.json(eval(response.SPE), expanded = False)

        except:
            st.write("failed to load chemspeed inventory")
            submitted = st.form_submit_button("re-fresh")

            if submitted:
                st.experimental_rerun()

    


            







