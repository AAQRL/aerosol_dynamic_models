# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\huangzhang\Desktop\Aerosol dynamics modeling\GUI\PyQt5 GUI\ADM\Dis_Dis-Sec\DisSec_Dialog_v1.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DisSec_Dialog(object):
    
    #HZ 2019  
    def __init__(self, InputTemp, InputPressure, InputSatPressure, \
                InputMW, InputDen, InputSurfTension):
        
        global Temp, Pressure, SatPressure, Den, MW, SurfT
        
        Temp = InputTemp
        Pressure = InputPressure
        SatPressure = InputSatPressure
        Den = InputDen
        MW = InputMW
        SurfT = InputSurfTension
        print("MW  = " + MW)
        print("Den  = " + Den)
                
    #HZ 2019   
    
    def setupUi(self, DisSec_Dialog):
        DisSec_Dialog.setObjectName("DisSec_Dialog")
        DisSec_Dialog.resize(551, 466)
        self.COMSOLHere = QtWidgets.QLabel(DisSec_Dialog)
        self.COMSOLHere.setGeometry(QtCore.QRect(10, 130, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.COMSOLHere.setFont(font)
        self.COMSOLHere.setObjectName("COMSOLHere")
        self.InputFinalSimuTime = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputFinalSimuTime.setGeometry(QtCore.QRect(440, 20, 104, 21))
        self.InputFinalSimuTime.setObjectName("InputFinalSimuTime")
        self.InputCond = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputCond.setGeometry(QtCore.QRect(200, 50, 104, 21))
        self.InputCond.setObjectName("InputCond")
        self.InputCOMSOL = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputCOMSOL.setGeometry(QtCore.QRect(200, 140, 104, 21))
        self.InputCOMSOL.setObjectName("InputCOMSOL")
        self.InputNucl = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputNucl.setGeometry(QtCore.QRect(200, 110, 104, 21))
        self.InputNucl.setObjectName("InputNucl")
        self.InputInitNumbConc = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputInitNumbConc.setGeometry(QtCore.QRect(200, 170, 104, 21))
        self.InputInitNumbConc.setObjectName("InputInitNumbConc")
        self.IniNumbConcHere = QtWidgets.QLabel(DisSec_Dialog)
        self.IniNumbConcHere.setGeometry(QtCore.QRect(10, 160, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.IniNumbConcHere.setFont(font)
        self.IniNumbConcHere.setObjectName("IniNumbConcHere")
        self.FinalSimuTimeHere = QtWidgets.QLabel(DisSec_Dialog)
        self.FinalSimuTimeHere.setGeometry(QtCore.QRect(310, 10, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FinalSimuTimeHere.setFont(font)
        self.FinalSimuTimeHere.setObjectName("FinalSimuTimeHere")
        self.CondHere = QtWidgets.QLabel(DisSec_Dialog)
        self.CondHere.setGeometry(QtCore.QRect(10, 40, 191, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CondHere.setFont(font)
        self.CondHere.setObjectName("CondHere")
        self.DisNumbHere = QtWidgets.QLabel(DisSec_Dialog)
        self.DisNumbHere.setGeometry(QtCore.QRect(10, 10, 111, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DisNumbHere.setFont(font)
        self.DisNumbHere.setObjectName("DisNumbHere")
        self.InputDisNumb = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputDisNumb.setGeometry(QtCore.QRect(200, 20, 104, 21))
        self.InputDisNumb.setObjectName("InputDisNumb")
        self.InputCoag = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputCoag.setGeometry(QtCore.QRect(200, 80, 104, 21))
        self.InputCoag.setObjectName("InputCoag")
        self.CoagHere = QtWidgets.QLabel(DisSec_Dialog)
        self.CoagHere.setGeometry(QtCore.QRect(10, 70, 181, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CoagHere.setFont(font)
        self.CoagHere.setObjectName("CoagHere")
        self.NuclHere = QtWidgets.QLabel(DisSec_Dialog)
        self.NuclHere.setGeometry(QtCore.QRect(10, 100, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NuclHere.setFont(font)
        self.NuclHere.setObjectName("NuclHere")
        self.Run = QtWidgets.QPushButton(DisSec_Dialog)
        self.Run.setGeometry(QtCore.QRect(210, 380, 75, 23))
        self.Run.setObjectName("Run")
        self.FinalSimuTimeHere_2 = QtWidgets.QLabel(DisSec_Dialog)
        self.FinalSimuTimeHere_2.setGeometry(QtCore.QRect(10, 200, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FinalSimuTimeHere_2.setFont(font)
        self.FinalSimuTimeHere_2.setObjectName("FinalSimuTimeHere_2")
        self.InputSizeFac = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputSizeFac.setGeometry(QtCore.QRect(200, 200, 104, 21))
        self.InputSizeFac.setObjectName("InputSizeFac")
        self.dmaxHere = QtWidgets.QLabel(DisSec_Dialog)
        self.dmaxHere.setGeometry(QtCore.QRect(10, 230, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.dmaxHere.setFont(font)
        self.dmaxHere.setObjectName("dmaxHere")
        self.Inputdmax = QtWidgets.QTextEdit(DisSec_Dialog)
        self.Inputdmax.setGeometry(QtCore.QRect(200, 240, 104, 21))
        self.Inputdmax.setObjectName("Inputdmax")
        self.TimeStage1Here = QtWidgets.QLabel(DisSec_Dialog)
        self.TimeStage1Here.setGeometry(QtCore.QRect(310, 50, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.TimeStage1Here.setFont(font)
        self.TimeStage1Here.setObjectName("TimeStage1Here")
        self.InputTimeStage1 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputTimeStage1.setGeometry(QtCore.QRect(440, 60, 104, 21))
        self.InputTimeStage1.setObjectName("InputTimeStage1")
        self.TimeStage2Here = QtWidgets.QLabel(DisSec_Dialog)
        self.TimeStage2Here.setGeometry(QtCore.QRect(310, 130, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.TimeStage2Here.setFont(font)
        self.TimeStage2Here.setObjectName("TimeStage2Here")
        self.NumbTimeStage1Here = QtWidgets.QLabel(DisSec_Dialog)
        self.NumbTimeStage1Here.setGeometry(QtCore.QRect(310, 90, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NumbTimeStage1Here.setFont(font)
        self.NumbTimeStage1Here.setObjectName("NumbTimeStage1Here")
        self.InputNumbTimeStage1 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputNumbTimeStage1.setGeometry(QtCore.QRect(440, 100, 104, 21))
        self.InputNumbTimeStage1.setObjectName("InputNumbTimeStage1")
        self.InputTimeStage2 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputTimeStage2.setGeometry(QtCore.QRect(440, 140, 104, 21))
        self.InputTimeStage2.setObjectName("InputTimeStage2")
        self.TimeStage3Here = QtWidgets.QLabel(DisSec_Dialog)
        self.TimeStage3Here.setGeometry(QtCore.QRect(310, 210, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.TimeStage3Here.setFont(font)
        self.TimeStage3Here.setObjectName("TimeStage3Here")
        self.NumbTimeStage2Here = QtWidgets.QLabel(DisSec_Dialog)
        self.NumbTimeStage2Here.setGeometry(QtCore.QRect(310, 170, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NumbTimeStage2Here.setFont(font)
        self.NumbTimeStage2Here.setObjectName("NumbTimeStage2Here")
        self.InputNumbTimeStage2 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputNumbTimeStage2.setGeometry(QtCore.QRect(440, 180, 104, 21))
        self.InputNumbTimeStage2.setObjectName("InputNumbTimeStage2")
        self.InputTimeStage3 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputTimeStage3.setGeometry(QtCore.QRect(440, 220, 104, 21))
        self.InputTimeStage3.setObjectName("InputTimeStage3")
        self.NumbTimeStage3 = QtWidgets.QLabel(DisSec_Dialog)
        self.NumbTimeStage3.setGeometry(QtCore.QRect(310, 250, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NumbTimeStage3.setFont(font)
        self.NumbTimeStage3.setObjectName("NumbTimeStage3")
        self.InputNumbTimeStage3 = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputNumbTimeStage3.setGeometry(QtCore.QRect(440, 260, 104, 21))
        self.InputNumbTimeStage3.setObjectName("InputNumbTimeStage3")
        self.NumbOutputHere = QtWidgets.QLabel(DisSec_Dialog)
        self.NumbOutputHere.setGeometry(QtCore.QRect(310, 290, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NumbOutputHere.setFont(font)
        self.NumbOutputHere.setObjectName("NumbOutputHere")
        self.InputNumbOutput = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputNumbOutput.setGeometry(QtCore.QRect(440, 300, 104, 21))
        self.InputNumbOutput.setObjectName("InputNumbOutput")
        self.EtaHere = QtWidgets.QLabel(DisSec_Dialog)
        self.EtaHere.setGeometry(QtCore.QRect(10, 270, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.EtaHere.setFont(font)
        self.EtaHere.setObjectName("EtaHere")
        self.InputEta = QtWidgets.QTextEdit(DisSec_Dialog)
        self.InputEta.setGeometry(QtCore.QRect(200, 280, 104, 21))
        self.InputEta.setObjectName("InputEta")

        self.retranslateUi(DisSec_Dialog)
        QtCore.QMetaObject.connectSlotsByName(DisSec_Dialog)
        
        #HZ 2019
        self.Run.clicked.connect(self.get_dissec_inputfile)
        #HZ 2019

    def retranslateUi(self, DisSec_Dialog):
        _translate = QtCore.QCoreApplication.translate
        DisSec_Dialog.setWindowTitle(_translate("DisSec_Dialog", "Dialog"))
        self.COMSOLHere.setText(_translate("DisSec_Dialog", "COMSOL Flag (0 or 1)"))
        self.InputFinalSimuTime.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCond.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCOMSOL.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputNucl.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.IniNumbConcHere.setText(_translate("DisSec_Dialog", "Inital Numb Conc. (#/m^3)"))
        self.FinalSimuTimeHere.setText(_translate("DisSec_Dialog", "Final Simu. Time (s)"))
        self.CondHere.setText(_translate("DisSec_Dialog", "Index for Condensation (0 or 1)"))
        self.DisNumbHere.setText(_translate("DisSec_Dialog", "Discrete Number"))
        self.CoagHere.setText(_translate("DisSec_Dialog", "Index for Coagulation (0 or 1)"))
        self.NuclHere.setText(_translate("DisSec_Dialog", "Nucleation Flag (0,1 or 2)"))
        self.Run.setText(_translate("DisSec_Dialog", "Run"))
        self.FinalSimuTimeHere_2.setText(_translate("DisSec_Dialog", "Size Factor"))
        self.InputSizeFac.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.dmaxHere.setText(_translate("DisSec_Dialog", "dmax (nm)"))
        self.Inputdmax.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.TimeStage1Here.setText(_translate("DisSec_Dialog", "Time for Stage 1 (S)"))
        self.InputTimeStage1.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.TimeStage2Here.setText(_translate("DisSec_Dialog", "Time for stage 2 (S)"))
        self.NumbTimeStage1Here.setText(_translate("DisSec_Dialog", "Numb. Time Stage 1"))
        self.InputNumbTimeStage1.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputTimeStage2.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.TimeStage3Here.setText(_translate("DisSec_Dialog", "Time for stage 3 (S)"))
        self.NumbTimeStage2Here.setText(_translate("DisSec_Dialog", "Numb. Time Stage 2"))
        self.InputNumbTimeStage2.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputTimeStage3.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.NumbTimeStage3.setText(_translate("DisSec_Dialog", "Numb. Time Stage 3"))
        self.InputNumbTimeStage3.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.NumbOutputHere.setText(_translate("DisSec_Dialog", "Number of Output"))
        self.InputNumbOutput.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.EtaHere.setText(_translate("DisSec_Dialog", "Eta"))
        self.InputEta.setHtml(_translate("DisSec_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
    
    #HZ 2019
    def get_dissec_inputfile(self):
        #Temp = self.InputTemp.toPlainText()
        #MW = self.InputMW.toPlainText()
        #Den = self.InputDen.toPlainText()
        #Temp = InputTemp
        #Den = InputDen
        #MW = InputMW
        global Temp, Pressure, SatPressure, Den, MW, SurfT
        print("Temp  = " + Temp)
        print("Pressure  = " + Pressure)
        print("SatPressure  = " + SatPressure)
        print("MW  = " + MW)
        print("Den  = " + Den)
        print("SurfT  = " + SurfT)
        
        #FD = self.InputFirstDiam.toPlainText() #first discrete diameter
        InitNumbConc = self.InputInitNumbConc.toPlainText()
        DisNumb = self.InputDisNumb.toPlainText()
        CoagFlag = self.InputCoag.toPlainText()
        CondFlag = self.InputCond.toPlainText()
        NuclFlag = self.InputNucl.toPlainText()
        SizeFac = self.InputSizeFac.toPlainText()
        dmax = self.Inputdmax.toPlainText()
        COMSOL = self.InputCOMSOL.toPlainText()
        #TotT = self.InputFinalSimuTIme.toPlainText()
        Time1 = self.InputTimeStage1.toPlainText()
        Time2 = self.InputTimeStage2.toPlainText()
        Time3 = self.InputTimeStage3.toPlainText()
        NumbTime1 = self.InputNumbTimeStage1.toPlainText()
        NumbTime2 = self.InputNumbTimeStage2.toPlainText()
        NumbTime3 = self.InputNumbTimeStage3.toPlainText()
        NumbOutput = self.InputNumbOutput.toPlainText()
        Eta = self.InputEta.toPlainText()
        #OutputonScreen = Temp + '\n' + MW + '\n' + Den + '\n' + FD + '\n' + InitNumbConc + '\n' \
        #                 + DisNumb + '\n' + CoagFlag + '\n' + TotT + '\n' + Dt
        
        #self.OutputScreen.setText(OutputonScreen)
        
        outfile = open('SYSOPH.DAT', 'w')
        
        MW1 = float(MW)*1000
        Den1 = float(Den)/1000
        SatPressure1 = float(SatPressure)*9.86923e-6
        SurfT1 = float(SurfT)*1000
        dmax1 = float(dmax)*10
        #print('MW1', MW1)        
        
        output_dissec = DisNumb + '\n' \
                 + str(0) + '\n' \
                 + str(MW1) + '\n' \
                 + str(Den1) + '\n' \
                 + str(0) + '\n' \
                 + str(0) + '\n' \
                 + InitNumbConc + '\n' \
                 + str(SatPressure1) + '\n' \
                 + str(SurfT1) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + ' ' + str(0) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + '\n' \
                 + str(0) + '\n' \
                 + str(143.7) + '\n' \
                 + str(5.7) + '\n' \
                 + str(0) + '\n' \
                 + str(0) + '\n' \
                 + str(0) + '\n' \
                 + str(0) + '\n' \
                 + str(307.0) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + ' ' + str(0) + '\n' \
                 + str(0) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + '\n' \
                 + str(1) + ' ' + str(0) + '\n' \
                 + SizeFac+ '\n' \
                 + str(dmax1) + '\n' \
                 + Temp + '\n' \
                 + COMSOL + '\n' \
                 + str(0) + '\n' \
                 + str(3) + '\n' \
                 + Time1 + '\n' \
                 + NumbTime1 + '\n' \
                 + Time2 + '\n' \
                 + NumbTime2 + '\n' \
                 + Time3 + '\n' \
                 + NumbTime3 + '\n' \
                 + NumbOutput + '\n' \
                 + NumbOutput + '\n' \
                 + str(30) + '\n' \
                 + str(4) + '\n' \
                 + Eta + '\n' \
                 + str(1) + '\n' \
                 + NuclFlag + '\n' \
                 + CondFlag + '\n' \
                 + CoagFlag + '\n' \
                 + str(0) + '\n' \
                 + 'concf.dat' + '\n' \
                 + 'output.dat'
                 
                 
        #CoagFlag + '\n' + DisNumb + '\n' + MW + ' ' + Den + ' ' + FD + \
        #         '\n' + InitNumbConc + '\n' + TotT + '\n' + Dt + '\n' + str(1e6) 
        
        outfile.write(output_dissec)
        
    #HZ 2019

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DisSec_Dialog = QtWidgets.QDialog()
    ui = Ui_DisSec_Dialog()
    ui.setupUi(DisSec_Dialog)
    DisSec_Dialog.show()
    sys.exit(app.exec_())

