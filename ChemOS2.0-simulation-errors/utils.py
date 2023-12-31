from sila2.client import SilaClient
from sila2.framework import SilaError
import time
import psycopg2 as psycopg
from sila2.client import SilaClient
import time
import json
from database import *
from sqlalchemy import (
    create_engine,
)
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import os
from pathlib import Path
import threading
from dblogin import get_login


dbname, user, dbpassword = get_login()
engine = create_engine(f"postgresql://{user}:{dbpassword}@localhost:5432/{dbname}")
database_connect = f"postgresql://{user}:{dbpassword}@localhost:5432/{dbname}"


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




def run_hplc_chemspeed(injection, ip, port, cert, id_dict):

    print("here")
    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{injection["name"]}',
                '{2}',
                '{json.dumps(injection)}',
                '{"test hplc job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    id_dict["id"] = JOB_ID

    def HPLC_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Data == "output_data":
            print("hplc data received")
            data = SubmitJob_response.Payload

            try:
                log = HPLCdata(
                name = "dummy hplc run",
                timestamp = jobtimestamp,
                code = data,
                job_id = JOB_ID
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Data)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{2}',
                                '{JOB_ID}',
                                '{SubmitJob_response.Data}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    #connect to hplcms
    client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.HPLCMSsimulator.SubmitJobChemspeed(str(injection))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(HPLC_update_SQLalchemy)

    while instance.done != True:
        time.sleep(1)

    return 

def hplc_from_autosampler(injection, ip, port, cert):

    print("here")
    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{injection["name"]}',
                '{2}',
                '{json.dumps(injection)}',
                '{"test hplc job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    def HPLC_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Data == "output_data":
            data = SubmitJob_response.Payload

            try:
                log = HPLCdata(
                name = injection["name"],
                timestamp = jobtimestamp,
                code = data,
                job_id = JOB_ID
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Data)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{2}',
                                '{JOB_ID}',
                                '{SubmitJob_response.Data}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    #connect to hplcms
    client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.HPLCMSsimulator.SubmitJobAutosampler(str(injection))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(HPLC_update_SQLalchemy)

    while instance.done != True:
        time.sleep(1)

    return JOB_ID

def hplc_from_chemspeed(injection, chemspeed_position, ip, port, cert, id_dict):

    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{injection["name"]}',
                '{2}',
                '{json.dumps(injection)}',
                '{"test hplc job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    id_dict["id"] = JOB_ID

    def HPLC_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Data == "output_data":
            data = SubmitJob_response.Payload

            try:
                log = HPLCdata(
                name = injection["name"],
                timestamp = jobtimestamp,
                code = data,
                job_id = JOB_ID
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Data)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{2}',
                                '{JOB_ID}',
                                '{SubmitJob_response.Data}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    #connect to hplcms
    if cert is None:
        client = SilaClient(ip, port, insecure=True) 
    else:
        client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.HPLCMSsimulator.SubmitJobChemspeed(str(injection))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(HPLC_update_SQLalchemy)


    while True:
        if HPLC_CERT is None:
            client2 = SilaClient(HPLC_IP,HPLC_PORT, insecure=True)
        else:
            client2 = SilaClient(HPLC_IP,HPLC_PORT, root_certs=open(HPLC_CERT, "rb").read())
        status = client2.HPLCMSsimulator.Status().Termination
        if status != "ready":
            print("waiting for ready signal from hplc")
            time.sleep(3)
        else:
            break
    print("injecting from chemspeed...")
    chemspeed_thread = threading.Thread(target=run_chemspeed_inject, args=(
        chemspeed_position,
        CHEMSPEED_IP, 
        CHEMSPEED_PORT, 
        CHEMSPEED_CERT
    ))
    chemspeed_thread.start()
    while chemspeed_thread.is_alive():
        time.sleep(1)
    while instance.done != True:
        time.sleep(1)
    if instance.status.value == 3:
        raise SilaError("HPLC has crashed")
    
    print("hplc thread done")

    return instance.get_responses().Termination


def run_chemspeed_filter(position, ip, port, cert):
    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, descriptions, timestamp, location) 
            VALUES (
                '{"filtration"}',
                '{1}',
                '{"filtration on chemspeed"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id



    def chemspeed_filter_update_SQLalchemy(filter_response) -> None:

        if filter_response.Status == "silasocket_operations":

            print("operations received")
            data = eval(filter_response.Operations)

            print(data)

            for operation in data:
                log = ChemspeedOperation(
                device_id =1,
                type_id = session.query(ChemspeedOperationType).filter_by(name=data[operation]["task"]).first().id,
                job_id = JOB_ID,
                code = json.dumps(data[operation]["parameters"]),
                timestamp = datetime.now(),
                location = "uoft"
                )       
                session.add(log) 
                session.commit()
            
        else:
            print(filter_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{1}',
                                '{JOB_ID}',
                                '{filter_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")

    if cert is None:
        client = SilaClient(ip,port, insecure=True)
    else:
        client = SilaClient(ip,port, root_certs=open(cert, "rb").read())
    
    instance = client.ChemSpeedOperator.Filter(position)

    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(chemspeed_filter_update_SQLalchemy)
    
    while instance.done != True:
        time.sleep(1)
    return instance

def run_chemspeed_inject(position, ip, port, cert):
    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, descriptions, timestamp, location) 
            VALUES (
                '{"Injection"}',
                '{1}',
                '{"Injection to chemspeed"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    def chemspeed_inject_update_SQLalchemy(inject_response) -> None:

        if inject_response.Status == "silasocket_operations":
            data = json.loads(inject_response.Operations)

            for operation in data:
                log = ChemspeedOperation(
                device_id =1,
                type_id = session.query(ChemspeedOperationType).filter_by(name=data[operation]["task"]).first().id,
                job_id = JOB_ID,
                code = json.dumps(data[operation]["parameters"]),
                timestamp = datetime.now(),
                location = "uoft"
                )       
                session.add(log) 
                session.commit()
            
        else:
            print(inject_response.Status)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{1}',
                                '{JOB_ID}',
                                '{inject_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")
    if cert is None:
        client = SilaClient(ip,port, insecure=True)
    else:
        client = SilaClient(ip,port, root_certs=open(cert, "rb").read())
    
    instance = client.ChemSpeedOperator.Inject(position)

    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(chemspeed_inject_update_SQLalchemy)
    
    while instance.done != True:
        time.sleep(1)
    return instance


def run_optics_table(jobdict, ip, port, cert):

    jobname = jobdict["name"]
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    jobtimestamp = datetime.now()


    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, descriptions, timestamp, location) 
            VALUES (
                '{"Injection"}',
                '{3}',
                '{"Optics table job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    job_id = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    def Optics_table_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Status == "output_data":

            print("output data received from optics table")
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
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
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

    if cert is None:
        client = SilaClient(ip, port, insecure=True)
    else:
        client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
    instance = client.OpticsTableSimulator.SubmitJob(str(jobdict))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(Optics_table_update_SQLalchemy)

    while instance.done != True:
        time.sleep(1)
    if instance.status.value == 3:
        raise SilaError("HPLC has crashed")

    return instance.get_responses().Termination, job_id

def hplc_blank():
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    jobtimestamp = datetime.now()

    with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, descriptions, timestamp, location) 
            VALUES (
                'blank',
                '{2}',
                '{"test hplc job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    def HPLC_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Data == "output_data":
            data = SubmitJob_response.Payload

            try:
                log = HPLCdata(
                name = "blank_run",
                timestamp = jobtimestamp,
                code = data,
                job_id = JOB_ID
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Data)
            try:
                with psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{2}',
                                '{JOB_ID}',
                                '{SubmitJob_response.Data}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")


    #connect to hplcms

    client = SilaClient(HPLC_IP, HPLC_PORT, insecure=True)
    instance = client.HPLCMSsimulator.BlankRun()
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(HPLC_update_SQLalchemy)

    while instance.done != True:
        time.sleep(1)

    return instance.get_responses()


    


def fetch_hplc_results(job_id):
    root = Path("jobresults/hplc/")
    conn = psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}")
    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT code FROM hplcdata WHERE job_id={job_id}")
            data,  = cur.fetchall()[0]
        if not os.path.exists(f"{root}/{job_id}"):
            os.mkdir(f"{root}/{job_id}")
        if data is None:
            return False
            
        with open(f"{root}/{job_id}/{job_id}.7z", "wb") as f:
            f.write(data)

        confirmfile = Path(f"{root}/{job_id}/confirm")

        if confirmfile.is_dir():
            return True

        else:
            os.system(f"7za e {root}/{job_id}/{job_id}.7z -aos -o{root}/{job_id}")
            os.mkdir(f"{root}/{job_id}/confirm")   
        return True
    except:
       return False

def fetch_optics_results(job_id):
    root = Path("jobresults/optics/")
    conn = psycopg.connect(f"dbname={dbname} user={user} password={dbpassword}")
    
    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT code FROM opticsdata WHERE job_id={job_id}")
            data,  = cur.fetchall()[0]
        if not os.path.exists(f"{root}/{job_id}"):
            os.mkdir(f"{root}/{job_id}")  
        if data is None:
            return False  
        with open(f"{root}/{job_id}/{job_id}.7z", "wb") as f:
            f.write(data)
        confirmfile = Path(f"{root}/{job_id}/confirm")

        if confirmfile.is_dir():
            return True
        else:
            os.system(f"7za e {root}/{job_id}/{job_id}.7z -aos -o{root}/{job_id}")
            os.mkdir(f"{root}/{job_id}/confirm")   
        return True
    except:
       return False
    

