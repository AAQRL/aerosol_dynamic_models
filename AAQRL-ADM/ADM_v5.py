# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\huangzhang\Desktop\Aerosol dynamics modeling\GUI\PyQt5 GUI\ADM\ADM_v2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
#from Dis_Dialog import Ui_Dis_Dialog
from Dis_Dialog_v2 import Ui_Dis_Dialog
from DisSec_Dialog_v1 import Ui_DisSec_Dialog
from Moment_Dialog import Ui_Moment_Dialog
from Modal_Dialog import Ui_Modal_Dialog

class Ui_MainWindow(object):
    
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(536, 357)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.SelctDis = QtWidgets.QPushButton(self.centralwidget)
        self.SelctDis.setGeometry(QtCore.QRect(430, 90, 75, 23))
        self.SelctDis.setObjectName("SelctDis")
        self.InputTemp = QtWidgets.QTextEdit(self.centralwidget)
        self.InputTemp.setGeometry(QtCore.QRect(180, 100, 104, 21))
        self.InputTemp.setObjectName("InputTemp")
        self.TempHere = QtWidgets.QLabel(self.centralwidget)
        self.TempHere.setGeometry(QtCore.QRect(10, 100, 131, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.TempHere.setFont(font)
        self.TempHere.setObjectName("TempHere")
        self.MWHere = QtWidgets.QLabel(self.centralwidget)
        self.MWHere.setGeometry(QtCore.QRect(10, 190, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.MWHere.setFont(font)
        self.MWHere.setObjectName("MWHere")
        self.InputMW = QtWidgets.QTextEdit(self.centralwidget)
        self.InputMW.setGeometry(QtCore.QRect(180, 200, 104, 21))
        self.InputMW.setObjectName("InputMW")
        self.DenHere = QtWidgets.QLabel(self.centralwidget)
        self.DenHere.setGeometry(QtCore.QRect(10, 230, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DenHere.setFont(font)
        self.DenHere.setObjectName("DenHere")
        self.InputDen = QtWidgets.QTextEdit(self.centralwidget)
        self.InputDen.setGeometry(QtCore.QRect(180, 240, 104, 21))
        self.InputDen.setObjectName("InputDen")
        self.BascParaHere = QtWidgets.QLabel(self.centralwidget)
        self.BascParaHere.setGeometry(QtCore.QRect(80, 60, 171, 20))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.BascParaHere.setFont(font)
        self.BascParaHere.setObjectName("BascParaHere")
        self.SelectModels = QtWidgets.QLabel(self.centralwidget)
        self.SelectModels.setGeometry(QtCore.QRect(350, 60, 141, 20))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.SelectModels.setFont(font)
        self.SelectModels.setObjectName("SelectModels")
        self.DisHere = QtWidgets.QLabel(self.centralwidget)
        self.DisHere.setGeometry(QtCore.QRect(330, 80, 91, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DisHere.setFont(font)
        self.DisHere.setObjectName("DisHere")
        self.DisSecHere = QtWidgets.QLabel(self.centralwidget)
        self.DisSecHere.setGeometry(QtCore.QRect(330, 120, 91, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DisSecHere.setFont(font)
        self.DisSecHere.setObjectName("DisSecHere")
        self.MomentHere = QtWidgets.QLabel(self.centralwidget)
        self.MomentHere.setGeometry(QtCore.QRect(330, 160, 91, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.MomentHere.setFont(font)
        self.MomentHere.setObjectName("MomentHere")
        self.ModalHere = QtWidgets.QLabel(self.centralwidget)
        self.ModalHere.setGeometry(QtCore.QRect(330, 200, 91, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ModalHere.setFont(font)
        self.ModalHere.setObjectName("ModalHere")
        self.AAQRL_ADM = QtWidgets.QLabel(self.centralwidget)
        self.AAQRL_ADM.setGeometry(QtCore.QRect(200, 20, 141, 20))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.AAQRL_ADM.setFont(font)
        self.AAQRL_ADM.setObjectName("AAQRL_ADM")
        self.PressureHere = QtWidgets.QLabel(self.centralwidget)
        self.PressureHere.setGeometry(QtCore.QRect(10, 130, 131, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.PressureHere.setFont(font)
        self.PressureHere.setObjectName("PressureHere")
        self.InputPressure = QtWidgets.QTextEdit(self.centralwidget)
        self.InputPressure.setGeometry(QtCore.QRect(180, 130, 104, 21))
        self.InputPressure.setObjectName("InputPressure")
        self.SatPressureHere = QtWidgets.QLabel(self.centralwidget)
        self.SatPressureHere.setGeometry(QtCore.QRect(10, 160, 151, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.SatPressureHere.setFont(font)
        self.SatPressureHere.setObjectName("SatPressureHere")
        self.InputSatPressure = QtWidgets.QTextEdit(self.centralwidget)
        self.InputSatPressure.setGeometry(QtCore.QRect(180, 170, 104, 21))
        self.InputSatPressure.setObjectName("InputSatPressure")
        self.SurfTensionHere = QtWidgets.QLabel(self.centralwidget)
        self.SurfTensionHere.setGeometry(QtCore.QRect(10, 260, 201, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.SurfTensionHere.setFont(font)
        self.SurfTensionHere.setObjectName("SurfTensionHere")
        self.InputSurfTension = QtWidgets.QTextEdit(self.centralwidget)
        self.InputSurfTension.setGeometry(QtCore.QRect(180, 270, 104, 21))
        self.InputSurfTension.setObjectName("InputSurfTension")
        self.SelctDis_Sec = QtWidgets.QPushButton(self.centralwidget)
        self.SelctDis_Sec.setGeometry(QtCore.QRect(430, 130, 75, 23))
        self.SelctDis_Sec.setObjectName("SelctDis_Sec")
        self.SelctMoment = QtWidgets.QPushButton(self.centralwidget)
        self.SelctMoment.setGeometry(QtCore.QRect(430, 170, 75, 23))
        self.SelctMoment.setObjectName("SelctMoment")
        self.SelctModal = QtWidgets.QPushButton(self.centralwidget)
        self.SelctModal.setGeometry(QtCore.QRect(430, 210, 75, 23))
        self.SelctModal.setObjectName("SelctModal")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 536, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave = QtWidgets.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
        #HZ 2019  
        self.SelctDis.clicked.connect(self.open_Dis_dialog)
        self.SelctDis_Sec.clicked.connect(self.open_DisSec_dialog)
        self.SelctMoment.clicked.connect(self.open_Moment_dialog)
        self.SelctModal.clicked.connect(self.open_Modal_dialog)
        
        #HZ 2019  

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "AAQRL-ADM"))
        self.SelctDis.setText(_translate("MainWindow", "Run"))
        self.InputTemp.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.TempHere.setText(_translate("MainWindow", "Temperature (K)"))
        self.MWHere.setText(_translate("MainWindow", "Molecualr Weight (kg/mol)"))
        self.DenHere.setText(_translate("MainWindow", "Particle Density (kg/m^3)"))
        self.BascParaHere.setText(_translate("MainWindow", "Species/Gas Properties"))
        self.SelectModels.setText(_translate("MainWindow", "Select Models"))
        self.DisHere.setText(_translate("MainWindow", "Discrete Model"))
        self.DisSecHere.setText(_translate("MainWindow", "Dis-Sec Model"))
        self.MomentHere.setText(_translate("MainWindow", "Moment Model"))
        self.ModalHere.setText(_translate("MainWindow", "Modal Model"))
        self.AAQRL_ADM.setText(_translate("MainWindow", "AAQRL-ADM"))
        self.PressureHere.setText(_translate("MainWindow", "Pressure (Pa)"))
        self.InputPressure.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.SatPressureHere.setText(_translate("MainWindow", "Saturation Pressure (Pa)"))
        self.InputSatPressure.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.SurfTensionHere.setText(_translate("MainWindow", "Particle Surf. Tension (N/m)"))
        self.SelctDis_Sec.setText(_translate("MainWindow", "Run"))
        self.SelctMoment.setText(_translate("MainWindow", "Run"))
        self.SelctModal.setText(_translate("MainWindow", "Run"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionSave.setText(_translate("MainWindow", "Save"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
    
    #HZ 2019    
    def open_Dis_dialog(self):
        Dis_Dialog = QtWidgets.QDialog()
        self.InputTemp = self.InputTemp.toPlainText()
        self.InputPressure = self.InputPressure.toPlainText()
        self.InputSatPressure = self.InputSatPressure.toPlainText()
        self.InputMW = self.InputMW.toPlainText()
        self.InputDen = self.InputDen.toPlainText()
        self.InputSurfTension = self.InputSurfTension.toPlainText()
        #print('Temp  = ' + self.InputTemp)
        ui = Ui_Dis_Dialog(self.InputTemp, self.InputPressure, \
                           self.InputSatPressure, self.InputMW, self.InputDen, \
                           self.InputSurfTension)
        #ui = Ui_Dis_Dialog()
        ui.setupUi(Dis_Dialog)
        Dis_Dialog.show()
        Dis_Dialog.exec_()
        
    def open_DisSec_dialog(self):
        
        DisSec_Dialog = QtWidgets.QDialog()
        self.InputTemp = self.InputTemp.toPlainText()
        self.InputPressure = self.InputPressure.toPlainText()
        self.InputSatPressure = self.InputSatPressure.toPlainText()
        self.InputMW = self.InputMW.toPlainText()
        self.InputMW = self.InputMW
        self.InputDen = self.InputDen.toPlainText()
        self.InputSurfTension = self.InputSurfTension.toPlainText()
        ui = Ui_DisSec_Dialog(self.InputTemp, self.InputPressure, \
                           self.InputSatPressure, self.InputMW, self.InputDen, \
                           self.InputSurfTension)
        
        ui.setupUi(DisSec_Dialog)
        DisSec_Dialog.show()
        DisSec_Dialog.exec_()
        
    def open_Moment_dialog(self):
        
        Moment_Dialog = QtWidgets.QDialog()
        self.InputTemp = self.InputTemp.toPlainText()
        self.InputPressure = self.InputPressure.toPlainText()
        self.InputSatPressure = self.InputSatPressure.toPlainText()
        self.InputMW = self.InputMW.toPlainText()
        self.InputMW = self.InputMW
        self.InputDen = self.InputDen.toPlainText()
        self.InputSurfTension = self.InputSurfTension.toPlainText()
        ui = Ui_Moment_Dialog(self.InputTemp, self.InputPressure, \
                           self.InputSatPressure, self.InputMW, self.InputDen, \
                           self.InputSurfTension)
        
        ui.setupUi(Moment_Dialog)
        Moment_Dialog.show()
        Moment_Dialog.exec_()
    
    def open_Modal_dialog(self):
        
        Modal_Dialog = QtWidgets.QDialog()
        self.InputTemp = self.InputTemp.toPlainText()
        self.InputPressure = self.InputPressure.toPlainText()
        self.InputSatPressure = self.InputSatPressure.toPlainText()
        self.InputMW = self.InputMW.toPlainText()
        self.InputMW = self.InputMW
        self.InputDen = self.InputDen.toPlainText()
        self.InputSurfTension = self.InputSurfTension.toPlainText()
        ui = Ui_Modal_Dialog(self.InputTemp, self.InputPressure, \
                           self.InputSatPressure, self.InputMW, self.InputDen, \
                           self.InputSurfTension)
        ui.setupUi(Modal_Dialog)
        Modal_Dialog.show()
        Modal_Dialog.exec_()
        
    #HZ 2019  
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

