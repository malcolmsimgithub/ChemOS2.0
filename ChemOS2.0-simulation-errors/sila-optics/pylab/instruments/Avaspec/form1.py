# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'form1.ui'
#
# Created by: PyQt5 UI code generator 5.8.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(679, 462)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.OpenCommBtn = QtWidgets.QPushButton(self.centralwidget)
        self.OpenCommBtn.setGeometry(QtCore.QRect(30, 50, 141, 23))
        self.OpenCommBtn.setObjectName("OpenCommBtn")
        self.CloseCommBtn = QtWidgets.QPushButton(self.centralwidget)
        self.CloseCommBtn.setGeometry(QtCore.QRect(30, 80, 141, 23))
        self.CloseCommBtn.setObjectName("CloseCommBtn")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(30, 150, 141, 191))
        self.groupBox.setObjectName("groupBox")
        self.IntTimeEdt = QtWidgets.QLineEdit(self.groupBox)
        self.IntTimeEdt.setGeometry(QtCore.QRect(10, 50, 51, 20))
        self.IntTimeEdt.setObjectName("IntTimeEdt")
        self.NumAvgEdt = QtWidgets.QLineEdit(self.groupBox)
        self.NumAvgEdt.setGeometry(QtCore.QRect(10, 100, 51, 20))
        self.NumAvgEdt.setObjectName("NumAvgEdt")
        self.NumMeasEdt = QtWidgets.QLineEdit(self.groupBox)
        self.NumMeasEdt.setGeometry(QtCore.QRect(10, 150, 51, 20))
        self.NumMeasEdt.setObjectName("NumMeasEdt")
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setGeometry(QtCore.QRect(10, 30, 111, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setGeometry(QtCore.QRect(10, 80, 111, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setGeometry(QtCore.QRect(10, 130, 131, 16))
        self.label_3.setObjectName("label_3")
        self.StartMeasBtn = QtWidgets.QPushButton(self.centralwidget)
        self.StartMeasBtn.setGeometry(QtCore.QRect(30, 380, 141, 23))
        self.StartMeasBtn.setObjectName("StartMeasBtn")
        self.StopMeasBtn = QtWidgets.QPushButton(self.centralwidget)
        self.StopMeasBtn.setGeometry(QtCore.QRect(30, 410, 141, 23))
        self.StopMeasBtn.setObjectName("StopMeasBtn")
        self.plot = RenderArea(self.centralwidget)
        self.plot.setGeometry(QtCore.QRect(190, 50, 471, 391))
        self.plot.setObjectName("plot")
        self.VersionBtn = QtWidgets.QPushButton(self.centralwidget)
        self.VersionBtn.setGeometry(QtCore.QRect(30, 110, 141, 23))
        self.VersionBtn.setObjectName("VersionBtn")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        MainWindow.setTabOrder(self.OpenCommBtn, self.CloseCommBtn)
        MainWindow.setTabOrder(self.CloseCommBtn, self.VersionBtn)
        MainWindow.setTabOrder(self.VersionBtn, self.IntTimeEdt)
        MainWindow.setTabOrder(self.IntTimeEdt, self.NumAvgEdt)
        MainWindow.setTabOrder(self.NumAvgEdt, self.NumMeasEdt)
        MainWindow.setTabOrder(self.NumMeasEdt, self.StartMeasBtn)
        MainWindow.setTabOrder(self.StartMeasBtn, self.StopMeasBtn)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.OpenCommBtn.setText(_translate("MainWindow", "Open Communication"))
        self.CloseCommBtn.setText(_translate("MainWindow", "Close Communication"))
        self.groupBox.setTitle(_translate("MainWindow", "Measurement Parameters"))
        self.label.setText(_translate("MainWindow", "Integration Time [ms]"))
        self.label_2.setText(_translate("MainWindow", "Number of Averages"))
        self.label_3.setText(_translate("MainWindow", "Number of Measurements"))
        self.StartMeasBtn.setText(_translate("MainWindow", "Start Measurements"))
        self.StopMeasBtn.setText(_translate("MainWindow", "Stop Measurements"))
        self.VersionBtn.setText(_translate("MainWindow", "Show Version Info"))

from renderarea import RenderArea
