# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from avaspec import *
import globals

class RenderArea(QWidget):
    points = QPolygonF(4096)
    def __init__(self, parent=None):
        super(RenderArea, self).__init__(parent)
        self.pen = QPen()
        self.pen.setCosmetic(True)  # line width must be independent of scale
        self.brush = QBrush()
        self.setBackgroundRole(QPalette.Base)
        self.setAutoFillBackground(True)
        # self.setAttribute(QtCore.Qt.WA_OpaquePaintEvent, False) makes no difference
#        self.timer = QTimer()
#        self.timer.timeout.connect(self.updatescreen)
#        self.timer.start(100) # 100 msec
        return

    #def minimumSizeHint(self):
    #    return QSize(500, 400)

    #def sizeHint(self):
    #    return QSize(500, 400)

    def setPen(self, pen):
        self.pen = pen
        self.update()

    def setBrush(self, brush):
        self.brush = brush
        self.update() 

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setPen(self.pen)
        painter.setBrush(self.brush)
        self.points.clear()
        x = 0
        while (x < globals.pixels):   # 0 through 2047
           self.points.append(QPointF(float(x), float(65536.0 - globals.spectraldata[x])))
           #self.points.append(QPointF(float(x), float(globals.spectraldata[x])))
           x += 1
        painter.scale(self.width()/globals.pixels, self.height()/65536.0)
        painter.drawPolyline(self.points)
        painter.setPen(self.palette().dark().color()) 
        painter.setBrush(Qt.NoBrush)
        painter.drawRect(0, 0, globals.pixels - 1, 65535)
        return

#    @pyqtSlot()
#    def updatescreen(self):
#        self.update()
