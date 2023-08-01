from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    ForeignKey,
    DateTime,
    LargeBinary
)

from sqlalchemy.dialects.mysql import LONGTEXT
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSON
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.functions import  current_timestamp
from datetime import datetime
# create the engine and connect to the PostgreSQL database
engine = create_engine("postgresql://postgres:postgres@localhost:5432/ifogchem")

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


# class Task(Base):
#     __tablename__ = "task"
#     id = Column(Integer, primary_key=True)
#     name = Column(String(64), unique=False)
#     timestamp = Column(DateTime)
#     # relation to device_log
#     device_log = relationship("DeviceLog", backref="task")

#     # repr
#     def __repr__(self):
#         return "<TaskLog %r>" % self.name



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

class PotentiostatData(Base):
    __tablename__ = "potentiostatdata"
    id = Column(Integer, primary_key=True)
    timestamp = Column(DateTime)
    return_data = Column(LONGTEXT)
    job_id = Column(Integer, ForeignKey("job.id"))
    def __repr__(self):
        return "<TaskLog %r>" % self.name



Base.metadata.create_all(engine)

# read and write
Session = sessionmaker(bind=engine)


