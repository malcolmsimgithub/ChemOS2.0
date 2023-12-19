from sila2.client import SilaClient
import time
import psycopg2 as psycopg
from sila2.client import SilaClient
import time
import json
from database import *
from sqlalchemy import (
    create_engine,
)
from dblogin import get_login
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import streamlit as st

dbname, dbuser, dbpassword = get_login()

database_connect = f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}"

def run_atlas_calculation(olympus, config, job_id, ip, port):

    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(f"postgresql://{dbuser}:{dbuser}@localhost:5432/{dbname}")
    Session = sessionmaker(bind=engine)
    session = Session()

    def printupdate(update) -> None:
        print(update.Status)
        if update.Status =='reccomended params':
            try:            
                session.query(AtlasData).filter(AtlasData.job_id == job_id).update({'return_data': update.StringPayload})
                session.commit()
            except:
                print("failed to log data recommended params")
        elif update.Status == "campaign_file":
            try:
                with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""UPDATE atlasdata SET localfile='{update.StringPayload}' WHERE(job_id={job_id});""")
                    conn.commit() 
            except:
                print("failed to log data campaign file")
        try:
            with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                with conn.cursor() as cur:
                        cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                        VALUES (
                            '{3}',
                            '{job_id}',
                            '{update.Status}',
                            '{jobtimestamp}',
                            '{"uoft"}'    
                            );
                        """)
                conn.commit() 
        except:
            print("failed to log data routine")
    
    client = SilaClient(ip, port, insecure=True)
    instance=client.Atlas.Recommend(olympus, str(config))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(printupdate)

    return instance



def run_chemspeed_synthesis(jsondata, jobname, job_id, ip, port,cert):

    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    jobtimestamp = datetime.now()
    if cert is None:
        print("insecure mode")
        client = SilaClient(ip,port, insecure=True)
    else:
        client = SilaClient(ip,port, root_certs=open(cert, "rb").read())

    def synthesis_update_SQLalchemy(addbatch_response) -> None:

        if addbatch_response.Status == "silasocket_operations":
            data = json.loads(addbatch_response.Operations)

            for operation in data:
                log = ChemspeedOperation(
                device_id =1,
                type_id = session.query(ChemspeedOperationType).filter_by(name=data[operation]["task"]).first().id,
                job_id = job_id,
                code = json.dumps(data[operation]["parameters"]),
                timestamp = datetime.now(),
                location = "uoft"
                )       
                session.add(log) 
                session.commit()
            
        else:
            print(addbatch_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{1}',
                                '{job_id}',
                                '{addbatch_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")
        
    jobdata = str(jsondata)
    
    instance = client.ChemSpeedOperator.Addbatch(jobname, jobdata)
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(synthesis_update_SQLalchemy) 

    return instance

    

def run_hplc_job(injection, jobname, job_id, ip, port, cert):
    
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    def HPLC_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Status == "output_data":
            print("check")
            data = SubmitJob_response.Payload

            try:
                log = HPLCdata(
                name = "dummy hplc run",
                timestamp = datetime.now(),
                code = data,
                job_id = job_id
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{2}',
                                '{job_id}',
                                '{SubmitJob_response.Status}',
                                '{datetime.now()}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    #connect to hplcms
    client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.HPLCMSsimulator.SubmitJob(str(injection))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(HPLC_update_SQLalchemy)

    return instance

def run_chemspeed_inject(position, name, smiles, filter, job_id, ip, port, cert):
    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    def chemspeed_inject_update_SQLalchemy(addbatch_response) -> None:

        if addbatch_response.Status == "silasocket_operations":
            data = json.loads(addbatch_response.Operations)

            for operation in data:
                log = ChemspeedOperation(
                device_id =1,
                type_id = session.query(ChemspeedOperationType).filter_by(name=data[operation]["task"]).first().id,
                job_id = job_id,
                code = json.dumps(data[operation]["parameters"]),
                timestamp = datetime.now(),
                location = "uoft"
                )       
                session.add(log) 
                session.commit()
            
        else:
            print(addbatch_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{1}',
                                '{job_id}',
                                '{addbatch_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    print("starting characterization, sending job to chemspeed for injection")
    client = SilaClient(ip,port, root_certs=open(cert, "rb").read())
    instance = client.ChemSpeedOperator.AddCharacterization(position, name, 
    smiles, 
    filter)
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(chemspeed_inject_update_SQLalchemy)
    
    return instance


def run_optics_table(jobdict, job_id, ip, port, cert):

    jobname = jobdict["name"]
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    jobtimestamp = datetime.now()


    def Optics_table_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Status == "output_data":
            data = SubmitJob_response.Payload
            try:
                log = OpticsData(
                name = jobname,
                timestamp = jobtimestamp,
                code = data,
                job_id = job_id
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={dbuser} password={dbuser}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{3}',
                                '{job_id}',
                                '{SubmitJob_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")

    client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.OpticsTableSimulator.SubmitJob(str(jobdict))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(Optics_table_update_SQLalchemy)

    return instance

def the_machine_job(data, job_id, ip, port):

    jobtimestamp = datetime.now()

    def printupdate(Submit_Job_Response) -> None:
        print(Submit_Job_Response.Status)
        try:
            with psycopg.connect(f"dbname={dbuser} user={dbuser} password={dbuser}") as conn:
                with conn.cursor() as cur:
                        cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                        VALUES (
                            '{5}',
                            '{job_id}',
                            '{Submit_Job_Response.Status}',
                            '{jobtimestamp}',
                            '{"uoft"}'    
                            );
                        """)
                conn.commit() 
        except:
            print("failed to log data")

    client = SilaClient("10.22.1.30", 65002, insecure=True)
    #client = SilaClient("127.0.0.1", 65002, insecure=True)

    instance = client.TheMachine.Runjob(data, "test.py")
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(printupdate)

    while instance.done != True:
        time.sleep(1)