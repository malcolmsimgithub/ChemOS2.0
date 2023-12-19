from . import CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import requests, time

import re
import numpy as np
from numpy.fft import fft, fftshift, fftfreq
import matplotlib.pyplot as plt

CMD_SHIM_MAP = CmdNameMap([
    (1, 'quick'),
    (2, 'medium'),
    (3, 'full'), ])

def str2fids(data):
    data = data.replace('\r', '').replace('\n\n', '\n')
    lmap = lambda f, l: list(map(f, l))

    # define regex rules
    header_re = re.compile(r'##([^\n=]+)=([^\n]+)')
    table_re  = re.compile(r'##PAGE=N=(\d)\n##([^\n]+)\n([^A-Za-z#]+)')
    value_re  = re.compile('-?[\d\.]+')

    headers = {k:v for k, v in header_re.findall(data)}
    facs    = lmap(float, headers['FACTOR'].split(','))
    zeroing = float(headers['$ZEROING'])

    fids = 0.j
    for match in table_re.findall(data):
        values = lmap(float, value_re.findall(match[2]))
        values = np.delete(values, slice(None, None, 5))
        if 'R' in match[1]:
            fids += values * facs[1]
        elif 'I' in match[1]:
            fids += 1.j * values * facs[2]
        else:
            print('There maybe errors')

    fids = np.append(fids, np.zeros(int(zeroing*len(fids))))
    return fids, headers

def phase_correction(spec, pc0, pc1):
    return spec * np.exp(1.0j * np.pi / 180.0 * (pc0 + pc1 * fftshift(fftfreq(len(spec), 1))))

def fids2spec(fids, headers, manual=True):
    dt       = float(headers['DELTAX'])
    ofreq    = float(headers['.OBSERVE FREQUENCY']) # in MHz
    N        = len(fids)
    o1       = float(headers['$O1']) # in Hz = SPECTRALCENTER * ofreq

    spec  = np.array( fftshift(fft(fids)) )

    # phase correction 
    if manual:
        pc0, pc1 = map(float, headers['$PHASECORRECTION'].split(','))
        spec = phase_correction(spec, pc0, pc1)

    idx   = fftshift(fftfreq(N, 1))
    hz    = idx/dt + o1
    ppm   = hz / ofreq
    return spec, hz, ppm


class NMR60Pro(object):
    def __init__(self, url = 'http://nmr.matterlab.sandbox:5000'):
        self.base_url = url + '/interfaces/'
        self._attr_path = []
        self.raw_data = None

    # call methods
    def fetch_json(self, res):
        if res.status_code == 200:
            return res.json()
        else:
            return {}

    def rest_get(self, path):
        res = requests.get(self.base_url+path)
        return self.fetch_json(res)

    def rest_put(self, path, data):
        res = requests.put(self.base_url+path, json = data)
        return self.fetch_json(res)

    # pipe methods
    def __getattr__(self, attr):
        self._attr_path.append(attr)
        return self

    def __call__(self, data = {}, **kwargs):
        # create url call path
        path = "/".join(self._attr_path)
        del self._attr_path[:]
        return self.call(path, data, **kwargs)

    def call(self, path, data = {}, **kwargs):
        req_data = data.copy()
        req_data.update(kwargs)
        if len(req_data.keys()) > 0:
            return self.rest_put(path, req_data)
        else:
            return self.rest_get(path)

    def set(self, **kwargs):
        return self.iFlow.ExperimentSettings(**kwargs)

    def get(self, *args):
        settings_dict = self.iFlow.ExperimentSettings()
        if len(args) == 0:
            return settings_dict
        elif len(args) == 1:
            return settings_dict.get(args[0], None)
        else:
            return [settings_dict.get(arg, None) for arg in args]

    # higher level experiment methods
    def idle(self):
        return self.iFlow.ExperimentStatus()['ResultCode'] == 0

    def run(self):
        self.iFlow.RunExperiment()
        while not self.idle():
            time.sleep(0.2)
        self.raw_data = self.iFlow.ExperimentStatus()['JDX_FileContents_TD']

    @mapsetmethod(CMD_SHIM_MAP)
    def shim(self, method):
        "Use lower case: quick, medium, full"
        self.iFlow.Shim({'ShimmingMethod': method, 'SolventShimming': False})
        while self.iFlow.Shim()['ShimmingMethod'] != 0:
            time.sleep(0.2)
        return self.iStatus.SpectrometerStatus()['Resolution']['LineWidths']

    def read_raw(self):
        if self.raw_data:
            return self.raw_data

    def read_fids(self):
        if self.raw_data:
            return str2fids(self.raw_data)[0]

    def read_spec(self):
        if self.raw_data:
            fids, headers = str2fids(self.raw_data)
            spec, _, ppm = fids2spec(fids, headers)
            return spec, ppm


