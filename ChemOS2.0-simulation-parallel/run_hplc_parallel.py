import asyncio

from sila2.client import SilaClient
import time
import time
import json
from database import *
from sqlalchemy import (
    create_engine,
)
from sqlalchemy.orm import sessionmaker
import json
from utils import *
from csv import reader
import sys


hplc_injection = {
        'name' : f"test_injection", 
        'itype' :'characterization_2nd',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : "CCC"},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

id_dict = {}

async def run_hplcs(injection, position):
    lock = asyncio.Lock()

    HPLC_IP_2 = "127.0.0.1"
    HPLC_PORT_2 = 65081

    HPLC_IP_1 = "127.0.0.1"
    HPLC_PORT_1 = 65010

    jobs = await asyncio.gather(
        hplc_from_chemspeed_async(injection, position, HPLC_IP_1, HPLC_PORT_1, None, lock),
        hplc_from_chemspeed_async(injection, position, HPLC_IP_2, HPLC_PORT_2, None, lock),
    )

    return jobs

result = asyncio.run(run_hplcs(hplc_injection, "SPE_C:9"))

if "not detected" in result:
    print("compound not detected. rerunning job")
    result = asyncio.run(run_hplcs(hplc_injection, "SPE_C:9"))
else:
    print("Compound detected in both samples!")
