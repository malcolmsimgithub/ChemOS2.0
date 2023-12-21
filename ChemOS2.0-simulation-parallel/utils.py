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
import asyncio


dbname, user, dbpassword = get_login()
engine = create_engine(f"postgresql://{user}:{dbpassword}@localhost:5432/{dbname}")
database_connect = f"postgresql://{user}:{dbpassword}@localhost:5432/{dbname}"


CHEMSPEED_IP = "127.0.0.1"
CHEMSPEED_CERT = None
CHEMSPEED_PORT = 65001
HPLC_PORT = 65010
OPTICS_IP = "127.0.0.1"
OPTICS_CERT = None
OPTICS_PORT = 65070



async def hplc_from_chemspeed_async(injection, chemspeed_position, ip, port, cert, lock):

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


    
    async with lock:
        #connect to hplcms
        if cert is None:
            client = SilaClient(ip, port, insecure=True) 
        else:
            client = SilaClient(ip, port, root_certs=open(cert, "rb").read())
        instance = client.HPLCMSsimulator.SubmitJobChemspeed(str(injection))
        sub = instance.subscribe_to_intermediate_responses()
        sub.add_callback(HPLC_update_SQLalchemy)
        while True:
            if cert is None:
                client2 = SilaClient(ip, port, insecure=True)
            else:
                client2 = SilaClient(ip, port, root_certs=open(cert, "rb").read())
            status = client2.HPLCMSsimulator.Status().Termination
            if status != "ready":
                print("waiting for ready signal from hplc")
                await asyncio.sleep(3)
            else:
                break
        print("injecting from chemspeed...")
        chemspeed_inject = await asyncio.to_thread(run_chemspeed_inject,
            chemspeed_position,
            CHEMSPEED_IP, 
            CHEMSPEED_PORT, 
            CHEMSPEED_CERT
        )
    while instance.done != True:
        await asyncio.sleep(3)
    if instance.status.value == 3:
        raise SilaError("HPLC has crashed")
    
    print("hplc thread done")

    return instance.get_responses().Termination

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
