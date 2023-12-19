from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    ForeignKey,
    DateTime,
    LargeBinary
)

from dblogin import get_login
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy.dialects.mysql import LONGTEXT
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSON
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.functions import  current_timestamp
from datetime import datetime
# create the engine and connect to the PostgreSQL database

dbname, dbuser, dbpassword = get_login()

engine = create_engine(f"postgresql://{dbuser}:{dbpassword}@localhost:5432/{dbname}")
if not database_exists(engine.url):
    create_database(engine.url)
# # create a declarative base class
Base = declarative_base()


class Device(Base):
    __tablename__ = "device"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    type = Column(String(64))
    descriptions = Column(String(64))
    manifacturing = Column(String(64))
    timestamp = Column(DateTime)
    location = Column(String(64))
    #device = relationship("Devices", backref="device", lazy=True)
    device_information = relationship("DeviceInformation", backref="device", lazy=True)
    job = relationship("Job", backref="device", lazy=True)
    device_log = relationship("DeviceLog", backref="device", lazy=True)
    # repr
    def __repr__(self):
        return "<Device %r>" % self.name


class DeviceInformation(Base):
    __tablename__ = "deviceinformation"
    id = Column(Integer, primary_key=True)
    device_id = Column(Integer, ForeignKey("device.id"))
    name = Column(String(64), unique=False)
    details = Column(String(64))
    version_driver = Column(String(64))
    timestamp = Column(DateTime)
    location = Column(String(64))

    # repr
    def __repr__(self):
        return "<DeviceInformation %r>" % (self.name)


class DeviceLog(Base):
    __tablename__ = "devicelog"
    id = Column(Integer, primary_key=True)
    device_id = Column(Integer, ForeignKey("device.id"))
    job_id = Column(Integer, ForeignKey("job.id"))
    name = Column(String(64), unique=False)
    status = Column(String(2000))
    timestamp = Column(DateTime)
    location = Column(String(64))

    # repre
    def __repr__(self):
        return "<DeviceLog %r>" % self.name


class Job(Base):
    __tablename__ = "job"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    devices_id = Column(Integer, ForeignKey("device.id"))
    code = Column(JSON)
    descriptions = Column(String(4000))
    timestamp = Column(DateTime)
    location = Column(String(64))

    # repr
    def __repr__(self):
        return "<Job %r>" % self.name


class ChemspeedOperation(Base):
    __tablename__ = "chemspeed_operation"
    id = Column(Integer, primary_key=True)
    device_id = Column(Integer, ForeignKey("device.id"))
    type_id = Column(Integer, ForeignKey("chemspeed_operation_type.id"))
    job_id = Column(Integer, ForeignKey("job.id"))
    name = Column(String(64), unique=False)
    details = Column(String(64))
    code = Column(JSON)
    timestamp = Column(DateTime)
    location = Column(String(64))

    # repr
    def __repr__(self):
        return "<Operation %r>" % self.name


class ChemspeedOperationType(Base):
    __tablename__ = "chemspeed_operation_type"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    options = Column(String(64))
    timestamp = Column(DateTime)
    # relation to oprations table
    operations = relationship("ChemspeedOperation", backref="chemspeed_operation_type")

    # repr
    def __repr__(self):
        return "<OperationType %r>" % self.name

class HPLCdata(Base):
    __tablename__ = "hplcdata"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    timestamp = Column(DateTime)
    code = Column(LargeBinary)
    job_id = Column(Integer, ForeignKey("job.id"))
    # repr
    def __repr__(self):
        return "<TaskLog %r>" % self.name

class OpticsData(Base):
    __tablename__ = "opticsdata"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    timestamp = Column(DateTime)
    code = Column(LargeBinary)
    job_id = Column(Integer, ForeignKey("job.id"))
    # repr
    def __repr__(self):
        return "<TaskLog %r>" % self.name

class AtlasData(Base):
    __tablename__ = "atlasdata"
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=False)
    timestamp = Column(DateTime)
    settings = Column(String(6000), unique=False)
    localfile = Column(String(100), unique=False)
    searchspace = Column(String(6000), unique=False)
    config = Column(JSON)
    campaign = Column(LargeBinary)
    return_data = Column(String(6000), unique = False)
    name = Column(String(64), unique=False)
    job_id = Column(Integer, ForeignKey("job.id"))
    # repr
    def __repr__(self):
        return "<TaskLog %r>" % self.name

Base.metadata.create_all(engine)

# read and write
Session = sessionmaker(bind=engine)
session = Session()

new_device = Device(
    name = "chemspeed",
    type = "synthesis",
    manifacturing = "chemspeed",
    timestamp = datetime.now(),
    location = "uoft"
)
session.add(new_device)
session.commit()

new_device = Device(
    name = "hplc",
    type = "purification",
    manifacturing = "thermo",
    timestamp = datetime.now(),
    location = "uoft"
)
session.add(new_device)
session.commit()

new_device = Device(
    name = "optics table",
    type = "characterization",
    manifacturing = "verious",
    timestamp = datetime.now(),
    location = "uoft"
)

session.add(new_device)
session.commit()


new_device = Device(
    name = "optimizer",
    type = "bayesian optimizer",
    manifacturing = "Riley Hickman",
    timestamp = datetime.now(),
    location = "chemos 2.0 internal"
)

session.add(new_device)
session.commit()

new_device = Device(
    name = "themachine",
    type = "synthesis",
    manifacturing = "various",
    timestamp = datetime.now(),
    location = "uoft"
)
session.add(new_device)
session.commit()



new_operation_type = ChemspeedOperationType(

    name = "transfer_compound",
    timestamp = datetime.now(),
)

session.add(new_operation_type)
session.commit()

new_operation_type = ChemspeedOperationType(

    name = "schlenk_cycle",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "reflux",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "prime_pumps",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "set_drawer",
    timestamp = datetime.now(),

)
session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "communicate",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "request_confirmation",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "filter_collect",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()
new_operation_type = ChemspeedOperationType(

    name = "submit_hplc_job",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()

new_operation_type = ChemspeedOperationType(

    name = "inject_to_hplc",
    timestamp = datetime.now(),

)

session.add(new_operation_type)
session.commit()


