#!/usr/bin/env python3
import os
import platform
import sys
import time
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from avaspec import *
import globals
import form1

class MainWindow(QMainWindow, form1.Ui_MainWindow):
    newdata = pyqtSignal()
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.IntTimeEdt.setText("{:3.1f}".format(5.0))
        self.NumAvgEdt.setText("{0:d}".format(1))
        self.NumMeasEdt.setText("{0:d}".format(1))
        self.StartMeasBtn.setEnabled(False)
        self.VersionBtn.setEnabled(False)
#       self.OpenCommBtn.clicked.connect(self.on_OpenCommBtn_clicked)
#       do not use explicit connect together with the on_ notation, or you will get
#       two signals instead of one!
        self.newdata.connect(self.handle_newdata)
   
    @pyqtSlot()
#   if you leave out the @pyqtSlot() line, you will also get an extra signal!
#   so you might even get three!
    def on_OpenCommBtn_clicked(self):
        ret = AVS_Init(0)    
        # QMessageBox.information(self,"Info","AVS_Init returned:  {0:d}".format(ret))
        ret = AVS_GetNrOfDevices()
        # QMessageBox.information(self,"Info","AVS_GetNrOfDevices returned:  {0:d}".format(ret))
        req = 0
        mylist = AvsIdentityType * 1
        ret = AVS_GetList(75, req, mylist)
        serienummer = str(ret[1].SerialNumber.decode("utf-8"))
        QMessageBox.information(self,"Info","Found Serialnumber: " + serienummer)
        globals.dev_handle = AVS_Activate(ret[1])
        # QMessageBox.information(self,"Info","AVS_Activate returned:  {0:d}".format(globals.dev_handle))
        devcon = DeviceConfigType
        reqsize = 0
        ret = AVS_GetParameter(globals.dev_handle, 63484, reqsize, devcon)
        globals.pixels = ret[1].m_Detector_m_NrPixels
        ret = AVS_GetLambda(globals.dev_handle,globals.wavelength)
        x = 0
        while (x < globals.pixels): # 0 through 2047
            globals.wavelength[x] = ret[x]
            x += 1
        self.StartMeasBtn.setEnabled(True)
        self.VersionBtn.setEnabled(True)
        return

    @pyqtSlot()
    def on_CloseCommBtn_clicked(self):
        callbackclass.callback(self, 0, 0)
        return

    @pyqtSlot()
    def on_VersionBtn_clicked(self):
        FPGAver = bytes(VERSION_LEN)
        FWver = bytes(VERSION_LEN)
        DLLver = bytes(VERSION_LEN)
        ret = AVS_GetVersionInfo(globals.dev_handle, FPGAver, FWver, DLLver)
        FPGAver = ret[0]
        FWver = ret[1]
        DLLver = ret[2]
        QMessageBox.information(self,"Info","FPGA version: {FPGA} \nFirmware version: {FW} \nDLL version: {DLL}" \
                               .format(FPGA=FPGAver.value.decode('utf-8'), 
                                       FW=FWver.value.decode('utf-8'),  
                                       DLL=DLLver.value.decode('utf-8')))
        return

    @pyqtSlot()
    def on_StartMeasBtn_clicked(self):
        ret = AVS_UseHighResAdc(globals.dev_handle, True)
        measconfig = MeasConfigType
        measconfig.m_StartPixel = 0
        measconfig.m_StopPixel = globals.pixels - 1
        measconfig.m_IntegrationTime = float(self.IntTimeEdt.text())
        measconfig.m_IntegrationDelay = 0
        measconfig.m_NrAverages = int(self.NumAvgEdt.text())
        measconfig.m_CorDynDark_m_Enable = 0  # nesting of types does NOT work!!
        measconfig.m_CorDynDark_m_ForgetPercentage = 0
        measconfig.m_Smoothing_m_SmoothPix = 0
        measconfig.m_Smoothing_m_SmoothModel = 0
        measconfig.m_SaturationDetection = 0
        measconfig.m_Trigger_m_Mode = 0
        measconfig.m_Trigger_m_Source = 0
        measconfig.m_Trigger_m_SourceType = 0
        measconfig.m_Control_m_StrobeControl = 0
        measconfig.m_Control_m_LaserDelay = 0
        measconfig.m_Control_m_LaserWidth = 0
        measconfig.m_Control_m_LaserWaveLength = 0.0
        measconfig.m_Control_m_StoreToRam = 0
        ret = AVS_PrepareMeasure(globals.dev_handle, measconfig)
        nummeas = int(self.NumMeasEdt.text())
        
        # to use Windows messages, supply a window handle to send the messages to
        # ret = AVS_Measure(globals.dev_handle, int(self.winId()), nummeas)
        # single message sent from DLL, confirmed with Spy++

        # when using polling, just pass a 0 for the windows handle
        scans = 0
        while (scans < nummeas):
            ret = AVS_Measure(globals.dev_handle, 0, 1)
            dataready = False
            while (dataready == False):
                dataready = (AVS_PollScan(globals.dev_handle) == True)
                time.sleep(0.001)
            if dataready == True:
                scans = scans + 1
                self.newdata.emit()
        return

    @pyqtSlot()
    def on_StopMeasBtn_clicked(self):
        ret = AVS_StopMeasure(globals.dev_handle)
        return

    #def nativeEvent(self, eventType, message):
    #    msg = ctypes.wintypes.MSG.from_address(message.__int__())
    #    if eventType == "windows_generic_MSG":
    #        if msg.message == WM_MEAS_READY:
    #            # print("Message Received!")
    #            self.newdata.emit()
    #    return False, 0

    @pyqtSlot()
    def handle_newdata(self):
        timestamp = 0
        ret = AVS_GetScopeData(globals.dev_handle, timestamp, globals.spectraldata )
        timestamp = ret[0]
        x = 0
        while (x < globals.pixels): # 0 through 2047
            globals.spectraldata[x] = ret[1][x]
            x += 1
            # QMessageBox.information(self,"Info","Received data")
        self.plot.update()
        return

def main():
    app = QApplication(sys.argv)
    app.lastWindowClosed.connect(app.quit)
    app.setApplicationName("PyQt5 simple demo")
    form = MainWindow()
    form.show()
    app.exec_()

main()
