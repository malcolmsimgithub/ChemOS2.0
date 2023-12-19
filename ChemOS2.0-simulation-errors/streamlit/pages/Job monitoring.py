import streamlit as st
import psycopg2
import pandas.io.sql as sqlio
import pandas as pd

from dblogin import get_login
dbname, dbuser, dbpassword = get_login()


@st.cache_resource
def init_connection():
    return psycopg2.connect(f"dbname={dbname} user={dbuser} password={dbpassword}")

conn = init_connection()

if st.button("refresh"):
    st.experimental_rerun()


# Perform query.
# Uses st.cache_data to only rerun when the query changes or after 10 min.
@st.cache_data(ttl=1)
def get_device():
    with conn.cursor() as cur:
        cur.execute("SELECT name, type, id, manifacturing  from device;")
        return cur.fetchall()

@st.cache_data(ttl=1)
def get_jobs(device, date):
    if len(device) == 1:
        device, = device
        with conn.cursor() as cur:
            cur.execute(f"SELECT name, id FROM job WHERE devices_id={device} AND timestamp BETWEEN '{date}T00:00:00.00' AND '{date}T23:59:59.999'")
            return cur.fetchall()   
    else:
        with conn.cursor() as cur:
            cur.execute(f"SELECT name, id FROM job WHERE devices_id IN {device} AND timestamp BETWEEN '{date}T00:00:00.00' AND '{date}T23:59:59.999'")
            return cur.fetchall()

@st.cache_data(ttl=1)
def get_job_info(job_id):
    with conn.cursor() as cur:
        cur.execute(f'SELECT status, timestamp, location from devicelog WHERE job_id={job_id}')
        return cur.fetchall()


st.subheader("Devices of the lab")
rows = st.dataframe(pd.DataFrame(get_device(), columns=["name", "purpose", "id", "manufacturing"]))

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


st.subheader("get job info")
with st.form("Jobinfo"):

    job_id = st.number_input("job_id- must be an integer!!", min_value=0, value=0, step=1)
    submitted = st.form_submit_button("Submit")
    if submitted:

        if job_id == 0:
            st.error("must be an integer value greater than zero")
        else:
            st.text("found job and loaded device_log")
            data = st.dataframe(pd.DataFrame(get_job_info(job_id), columns=["info", "timestamp", "location"]))
            



