import streamlit as st


from utils_streamlit import run_atlas_calculation


from database import *
import psycopg2 as psycopg

import json
from database import *
from sqlalchemy import (
    create_engine,

)
import os
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import psycopg2
import streamlit as st
import pandas.io.sql as sqlio
import pandas as pd


from dblogin import get_login
dbname, dbuser, dbpassword = get_login()

engine = create_engine(f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}")
database_connect = f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}"

ATLAS_IP = "127.0.0.1"

@st.cache(allow_output_mutation=True)
def init_connection():
    return psycopg2.connect(f"dbname={dbname} user={dbuser} password={dbpassword}")

conn = init_connection()

tab1, tab2 = st.tabs(["query", "get results"])


with tab1:
    st.title("Atlas Bayesian optimizer")
    with st.form("AtlasJob"):

        st.subheader("Use the Atlas Planner")
        st.write("Upload a json file of the planner configuration, or fill out the form")

        jobname = st.text_input("name your job", max_chars=100)
        description = st.text_input("job description", max_chars= 1000)
        optimizer_kind = st.selectbox("select optimizer type", options=["gradient", "genetic"])
        optimizer_type = st.selectbox("select acquisition function type", options=["EI", "PI", "UCB"])
        goal = st.selectbox("goal of campaign", options=["maximize", "minimize"])

        uploader = st.file_uploader("upload your planner config file, or enter into form", type="json")
        uploader2 = st.file_uploader("upload your Olympus Campaign", type="pkl")

        placeholder = st.empty()

        with placeholder.container():
            submitted = st.form_submit_button("Submit")
        
        if submitted:
            if jobname == None or uploader2 == None:
                st.text("jobname is invalid")
            
            olympusdata = uploader2.getvalue()



        if submitted:

            config = {}
            config["name"] = jobname
            config['goal'] = goal
            config['num_init_design'] = 1
            config['has_descriptors'] = False
            config['acquisition_optimizer_kind'] = optimizer_kind
            config['acquisition_type'] = optimizer_type

            print(config['acquisition_optimizer_kind'])

            engine = create_engine(database_connect)
            Session = sessionmaker(bind=engine)
            session = Session()

            jobtimestamp = datetime.now()

            with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbpassword}") as conn:
                with conn.cursor() as cur:
                    cur.execute(f"""INSERT INTO job(name, devices_id,descriptions, timestamp, location) 
                    VALUES (
                        '{jobname}',
                        '{3}',
                        '{description}',
                        '{jobtimestamp}',
                        '{"uoft"}'    
                        );
                    """)
                    conn.commit() 
                
            
            JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

            log = AtlasData(
                name = config["name"],
                job_id = JOB_ID,
                campaign = olympusdata,
                config = json.dumps(config),
                timestamp = datetime.now()
                )       
            session.add(log) 
            session.commit()

            st.write(f"your job id is {JOB_ID}")

            instance = run_atlas_calculation(olympusdata, config, JOB_ID, ATLAS_IP, 65100)


            st.warning("Job submitted successfully- to view results: go to display data. to submit a new job, refresh the page or click on the button below")

            placeholder.empty()

            submitted = st.form_submit_button("re-fresh and submit new job")

            if submitted:
                st.experimental_rerun()



with tab2:
    st.write("check on results of planner")

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
    def get_atlas_results(job_id):

        try:
            with conn.cursor() as cur:
                cur.execute(f"SELECT return_data FROM opticsdata WHERE job_id={job_id}")
                data,  = cur.fetchall()[0]
            
          

            return data


        except:
            return False
        
    st.subheader("See Jobs submitted")
    with st.form("Jobs"):

        df = sqlio.read_sql_query("select id, name from device", conn)

        device = st.multiselect("select devices", options=df["name"])
        date = st.date_input("choose date of job")

        submitted = st.form_submit_button("Submit")
        if submitted:

            
            devices = tuple(device)

            if len(devices) <1:
                st.error("please select at least one device")
            
            else:


                with conn.cursor() as cur:
                    command = 'SELECT id from device WHERE name IN %s'
                    cur.execute(command, (devices,))
                    id = [r[0] for r in cur.fetchall()]
        
                data = st.dataframe(pd.DataFrame(get_jobs(tuple(id), date), columns=["name", "job id"]))
                
    with st.form("Jobinfo"):

        job_id = st.number_input("job_id- must be an integer!!", min_value=0, value=0, step=1)
        submitted = st.form_submit_button("Submit")
        if submitted:

            if job_id == 0:
                st.error("must be an integer value greater than zero")
            else:
                st.text("found job")
                st.text("atlas reccomended params are:")
                st.write(get_atlas_results(job_id))
            
        
        
    
