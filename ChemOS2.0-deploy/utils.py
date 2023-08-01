from sila2.client import SilaClient
import time
import psycopg2 as psycopg
from sila2.client import SilaClient
import time
import json
from connection import *
from sqlalchemy import (
    create_engine,
)
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import json
import os
from pathlib import Path
import threading

user = "malcolm"
database_connect = "postgresql://malcolm:malcolm@localhost:5432/ifogchem"
"10.22.1.20"
CHEMSPEED_IP = "10.22.1.14"
CHEMSPEED_CERT = None
HPLC_IP = "10.22.1.21"
HPLC_CERT = "certificates/hplc-deploy.crt"
OPTICS_IP = "10.22.1.21"
OPTICS_CERT = "certificates/optics-deploy.crt"
ATLAS_IP = "127.0.0.1"
POTENTIOSTAT_IP="10.22.10.103"



def run_atlas_calculation(olympus, config, ip, port):

    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(f"postgresql://{user}:{user}@localhost:5432/ifogchem")
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                    with conn.cursor() as cur:
                        cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
                        VALUES (
                            '{config["jobname"]}',
                            '{4}',
                            '{json.dumps(config)}',
                            '{"test hplc job"}',
                            '{jobtimestamp}',
                            '{"uoft"}'    
                            );
                        """)
                        conn.commit() 
    
    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    log = AtlasData(
                name = config["jobname"],
                job_id = JOB_ID,
                campaign = olympus,
                config = json.dumps(config),
                timestamp = datetime.now()
                )       
    session.add(log) 
    session.commit()

    def printupdate(update) -> None:
        print(update.Status)
        if update.Status =='reccomended params':
            try:            
                session.query(AtlasData).filter(AtlasData.job_id == JOB_ID).update({'return_data': update.StringPayload})
                session.commit()
            except:
                print("failed to log data recommended params")
        elif update.Status == "campaign_file":
            try:
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""UPDATE atlasdata SET localfile='{update.StringPayload}' WHERE(job_id={JOB_ID});""")
                    conn.commit() 
            except:
                print("failed to log data campaign file")
        try:
            with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                with conn.cursor() as cur:
                        cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                        VALUES (
                            '{4}',
                            '{JOB_ID}',
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

    data = None
    while data == None:
        try:
            data = session.query(AtlasData).filter_by(job_id = JOB_ID).first().return_data
        except:
            print("waiting to retrieve reccomended params...")
        time.sleep(1)

    reccomends = eval(data)

    return reccomends


def run_chemspeed_synthesis(jsondata, jobname, ip, port,cert):

    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    jobtimestamp = datetime.now()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{jobname}',
                '{1}',
                '{json.dumps(jsondata)}',
                '{"another synthesis of atlas reccomendations"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    print(JOB_ID)

    def synthesis_update_SQLalchemy(addbatch_response) -> None:

        if addbatch_response.Status == "silasocket_operations":

            print("received")            
            data = json.loads(addbatch_response.Operations)

            for operation in data:
                log = ChemspeedOperation(
                device_id =1,
                type_id = session.query(ChemspeedOperationType).filter_by(name=data[operation]["task"]).first().id,
                job_id = JOB_ID,
                code = json.dumps(data[operation]),
                timestamp = datetime.now(),
                location = "uoft"
                )       
                session.add(log) 
                session.commit()
            
        else:
            print(addbatch_response.Status)
            try:
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{1}',
                                '{JOB_ID}',
                                '{addbatch_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")
        
    jobdata = str(jsondata)

    if cert is None:
        client = SilaClient(ip,port, insecure=True)
    else:
        client = SilaClient(ip,port, root_certs=open(cert, "rb").read())
    
    instance = client.ChemSpeedOperator.Addbatch(jobname, jobdata)
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(synthesis_update_SQLalchemy) 

    # while instance.done != True:
    #     time.sleep(1)

    return instance

    

def run_hplc_chemspeed(injection, ip, port, cert, id_dict):

    print("here")
    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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
            print("check")
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
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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

    print("hplc thread done")

    return 


def hplc_from_autosampler(injection, ip, port, cert):

    print("here")
    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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
            print("check")
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
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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

    print("hplc thread done")

    return JOB_ID

def hplc_from_chemspeed(injection, ip, port, cert, id_dict):

    print("here")
    jobtimestamp = datetime.now()
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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
            print("check")
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
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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

    print("hplc thread done")

    return JOB_ID


def run_chemspeed_filter(position, ip, port, cert):
    global database_connect
    jobtimestamp = datetime.now()
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()


    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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


    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
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

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{f"{jobname}"}',
                '{3}',
                '{json.dumps(jobdict)}',
                '{"test optics table job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id


    def Optics_table_update_SQLalchemy(SubmitJob_response) -> None:
        if SubmitJob_response.Status == "output_data":
            data = SubmitJob_response.Payload
            try:
                log = OpticsData(
                name = jobname,
                timestamp = jobtimestamp,
                code = data,
                job_id = JOB_ID
                )       
                session.add(log) 
                session.commit()
            except:
                print("error")
            
        else:
            print(SubmitJob_response.Status)
            try:
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{3}',
                                '{JOB_ID}',
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

    return instance, JOB_ID


def run_potentiostat(jobdict, ip, port, cert):

    jobname = jobdict["name"]
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    jobtimestamp = datetime.now()

    with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
        with conn.cursor() as cur:
            cur.execute(f"""INSERT INTO job(name, devices_id, code, descriptions, timestamp, location) 
            VALUES (
                '{f"{jobname}"}',
                '{6}',
                '{json.dumps(jobdict)}',
                '{"potentiostat job"}',
                '{jobtimestamp}',
                '{"uoft"}'    
                );
            """)
            conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    print(f"JOB ID IS {JOB_ID}")

    def potentiostat_update_SQLalchemy(SubmitJob_response) -> None:
        print(SubmitJob_response.Status)
        if SubmitJob_response.Status == "output_data":
            print("data received")
            data = SubmitJob_response.Payload
            
            with open("test.csv", "w") as f:
                f.write(SubmitJob_response.Payload)
            log = PotentiostatData(
            timestamp = jobtimestamp,
            return_data = data,
            job_id = JOB_ID
            )       
            session.add(log) 
            session.commit()
            print("CSV saved")
                # except:
                #     print("error saving csv")
            
        else:
            try:
                with psycopg.connect(f"dbname=ifogchem user={user} password={user}") as conn:
                    with conn.cursor() as cur:
                            cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                            VALUES (
                                '{6}',
                                '{JOB_ID}',
                                '{SubmitJob_response.Status}',
                                '{jobtimestamp}',
                                '{"uoft"}'    
                                );
                            """)
                    conn.commit() 
            except:
                print("failed to log data")

    client = SilaClient(ip, port, insecure=True)
    #instance = client.PotenServer.PrepareSample(json.dumps(jobdict))
    instance = client.PotenServer.RunExp(json.dumps(jobdict))
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(potentiostat_update_SQLalchemy)

    while instance.done !=True:
        time.sleep(1)

    time.sleep(10)

    redox_pot = instance.get_responses().Termination

    return JOB_ID, redox_pot


def run_themachine(name, jobfilepath, pyscriptname):
    global database_connect
    engine = create_engine(database_connect)
    Session = sessionmaker(bind=engine)
    session = Session()

    jobtimestamp = datetime.now()

    with psycopg.connect("dbname=ifogchem user=malcolm password=malcolm") as conn:
            with conn.cursor() as cur:
                cur.execute(f"""INSERT INTO job(name, devices_id, descriptions, timestamp, location) 
                VALUES (
                    '{name}',
                    '{5}',
                    '{""}',
                    '{jobtimestamp}',
                    '{"uoft"}'    
                    );
                """)
                conn.commit() 

    JOB_ID = session.query(Job).filter_by(timestamp=jobtimestamp).first().id

    def printupdate(Submit_Job_Response) -> None:
        print(Submit_Job_Response.Status)
        try:
            time = datetime.now()
            with psycopg.connect("dbname=ifogchem user=malcolm password=malcolm") as conn:
                with conn.cursor() as cur:
                        cur.execute(f"""INSERT INTO devicelog(device_id, job_id, status, timestamp, location) 
                        VALUES (
                            '{5}',
                            '{JOB_ID}',
                            '{Submit_Job_Response.Status}',
                            '{time}',
                            '{"uoft"}'    
                            );
                        """)
                conn.commit() 
        except:
            print("failed to log data")

    client = SilaClient("10.22.1.30", 65002, insecure=True)

    with open(jobfilepath, "r") as f:
        data = f.read()

    instance = client.TheMachine.Runjob(data, pyscriptname)
    sub = instance.subscribe_to_intermediate_responses()
    sub.add_callback(printupdate)

    while instance.done != True:
        time.sleep(1)



def hplc_and_chemspeed_inject(hplc_job, chemspeed_position):

    resultsdict = {}

    print("creating hplc job...")

    hplc_thread = threading.Thread(target=hplc_from_chemspeed, args=(hplc_job, HPLC_IP, 65004, HPLC_CERT, resultsdict))

    hplc_thread.start()

    while True:
        client = SilaClient(HPLC_IP,65004, root_certs=open(HPLC_CERT, "rb").read())
        status = client.HPLCMSsimulator.Status().Termination
        if status != "ready":
            print("waiting for ready signal from hplc")
            time.sleep(3)
        else:
            break
    
    print(resultsdict)

    print("injecting from chemspeed...")
    chemspeed_thread = threading.Thread(target=run_chemspeed_inject, args=(
        chemspeed_position,
        CHEMSPEED_IP, 
        65002, 
        CHEMSPEED_CERT
    ))

    chemspeed_thread.start()

    return hplc_thread, chemspeed_thread, resultsdict


def fetch_hplc_results(job_id):
    conn = psycopg.connect(f"dbname=ifogchem user=malcolm password=malcolm")
    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT code FROM hplcdata WHERE job_id={job_id}")
            data,  = cur.fetchall()[0]
        if not os.path.exists(f"bigjobresults/hplc/{job_id}"):
            os.mkdir(f"bigjobresults/hplc/{job_id}")
        if data is None:
            return False
         
        with open(f"bigjobresults/hplc/{job_id}/{job_id}.7z", "wb") as f:
            f.write(data)

        confirmfile = Path(f"bigjobresults/hplc/{job_id}/confirm")

        if confirmfile.is_dir():
            return True

        else:
            os.system(f"7za e bigjobresults/hplc/{job_id}/{job_id}.7z -aos -obigjobresults/hplc/{job_id}")
            os.mkdir(f"bigjobresults/hplc/{job_id}/confirm")   
        return True
    except:
       return False

def fetch_optics_results(job_id):
    root = Path("bigjobresults/optics/")
    conn = psycopg.connect(f"dbname=ifogchem user=malcolm password=malcolm")
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
    

