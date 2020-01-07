# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\huangzhang\Desktop\Aerosol dynamics modeling\GUI\PyQt5 GUI\ADM\Dis_Dis-Sec_Moment_Modal\Moment_Dialog.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Moment_Dialog(object):
    
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
    
    def setupUi(self, Moment_Dialog):
        Moment_Dialog.setObjectName("Moment_Dialog")
        Moment_Dialog.resize(543, 466)
        self.ReacHere = QtWidgets.QLabel(Moment_Dialog)
        self.ReacHere.setGeometry(QtCore.QRect(110, 250, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ReacHere.setFont(font)
        self.ReacHere.setObjectName("ReacHere")
        self.InputFinalSimuTIme = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputFinalSimuTIme.setGeometry(QtCore.QRect(280, 300, 104, 21))
        self.InputFinalSimuTIme.setObjectName("InputFinalSimuTIme")
        self.InputCond = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputCond.setGeometry(QtCore.QRect(280, 190, 104, 21))
        self.InputCond.setObjectName("InputCond")
        self.InputReac = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputReac.setGeometry(QtCore.QRect(280, 260, 104, 21))
        self.InputReac.setObjectName("InputReac")
        self.InputNucl = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputNucl.setGeometry(QtCore.QRect(280, 230, 104, 21))
        self.InputNucl.setObjectName("InputNucl")
        self.InputInitNumbConc = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputInitNumbConc.setGeometry(QtCore.QRect(280, 110, 104, 21))
        self.InputInitNumbConc.setObjectName("InputInitNumbConc")
        self.InputReactRate = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputReactRate.setGeometry(QtCore.QRect(280, 330, 104, 21))
        self.InputReactRate.setObjectName("InputReactRate")
        self.IniNumbConcHere = QtWidgets.QLabel(Moment_Dialog)
        self.IniNumbConcHere.setGeometry(QtCore.QRect(110, 100, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.IniNumbConcHere.setFont(font)
        self.IniNumbConcHere.setObjectName("IniNumbConcHere")
        self.FinalSimuTimeHere = QtWidgets.QLabel(Moment_Dialog)
        self.FinalSimuTimeHere.setGeometry(QtCore.QRect(110, 290, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FinalSimuTimeHere.setFont(font)
        self.FinalSimuTimeHere.setObjectName("FinalSimuTimeHere")
        self.CondHere = QtWidgets.QLabel(Moment_Dialog)
        self.CondHere.setGeometry(QtCore.QRect(90, 180, 181, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CondHere.setFont(font)
        self.CondHere.setObjectName("CondHere")
        self.segmagHere = QtWidgets.QLabel(Moment_Dialog)
        self.segmagHere.setGeometry(QtCore.QRect(110, 60, 111, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.segmagHere.setFont(font)
        self.segmagHere.setObjectName("segmagHere")
        self.Inputsegmag = QtWidgets.QTextEdit(Moment_Dialog)
        self.Inputsegmag.setGeometry(QtCore.QRect(280, 70, 104, 21))
        self.Inputsegmag.setObjectName("Inputsegmag")
        self.dpgHere = QtWidgets.QLabel(Moment_Dialog)
        self.dpgHere.setGeometry(QtCore.QRect(110, 30, 101, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.dpgHere.setFont(font)
        self.dpgHere.setObjectName("dpgHere")
        self.InputCoag = QtWidgets.QTextEdit(Moment_Dialog)
        self.InputCoag.setGeometry(QtCore.QRect(280, 150, 104, 21))
        self.InputCoag.setObjectName("InputCoag")
        self.CoagHere = QtWidgets.QLabel(Moment_Dialog)
        self.CoagHere.setGeometry(QtCore.QRect(110, 140, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CoagHere.setFont(font)
        self.CoagHere.setObjectName("CoagHere")
        self.ReactRateHere = QtWidgets.QLabel(Moment_Dialog)
        self.ReactRateHere.setGeometry(QtCore.QRect(90, 320, 171, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ReactRateHere.setFont(font)
        self.ReactRateHere.setObjectName("ReactRateHere")
        self.NuclHere = QtWidgets.QLabel(Moment_Dialog)
        self.NuclHere.setGeometry(QtCore.QRect(110, 220, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NuclHere.setFont(font)
        self.NuclHere.setObjectName("NuclHere")
        self.Inputdpg = QtWidgets.QTextEdit(Moment_Dialog)
        self.Inputdpg.setGeometry(QtCore.QRect(280, 40, 104, 21))
        self.Inputdpg.setObjectName("Inputdpg")
        self.RunMoment = QtWidgets.QPushButton(Moment_Dialog)
        self.RunMoment.setGeometry(QtCore.QRect(210, 380, 75, 23))
        self.RunMoment.setObjectName("RunMoment")

        self.retranslateUi(Moment_Dialog)
        QtCore.QMetaObject.connectSlotsByName(Moment_Dialog)
        
        #HZ 2019
        self.RunMoment.clicked.connect(self.get_Moment_inputfile)
        #HZ 2019

    def retranslateUi(self, Moment_Dialog):
        _translate = QtCore.QCoreApplication.translate
        Moment_Dialog.setWindowTitle(_translate("Moment_Dialog", "Dialog"))
        self.ReacHere.setText(_translate("Moment_Dialog", "Reac Flag (0 or 1)"))
        self.InputFinalSimuTIme.setHtml(_translate("Moment_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCond.setHtml(_translate("Moment_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputReac.setHtml(_translate("Moment_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputNucl.setHtml(_translate("Moment_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.IniNumbConcHere.setText(_translate("Moment_Dialog", "Inital Numb Conc. (#/m^3)"))
        self.FinalSimuTimeHere.setText(_translate("Moment_Dialog", "Final Simu. Time (s)"))
        self.CondHere.setText(_translate("Moment_Dialog", "Condensation Flag (0, 1 or 2)"))
        self.segmagHere.setText(_translate("Moment_Dialog", "Sigmag"))
        self.dpgHere.setText(_translate("Moment_Dialog", "dpg (m)"))
        self.CoagHere.setText(_translate("Moment_Dialog", "Coagulation  Flag (0 or 1)"))
        self.ReactRateHere.setText(_translate("Moment_Dialog", "Dimensionless reaction rate"))
        self.NuclHere.setText(_translate("Moment_Dialog", "Nucleation Flag (0, 1 or 2)"))
        self.RunMoment.setText(_translate("Moment_Dialog", "Run"))


    def get_Moment_inputfile(self):
        
        global Temp, Pressure, SatPressure, Den, MW, SurfT
        print("Temp  = " + Temp)
        print("Pressure  = " + Pressure)
        print("SatPressure  = " + SatPressure)
        print("MW  = " + MW)
        print("Den  = " + Den)
        print("SurfT  = " + SurfT)
              
        dpg = self.Inputdpg.toPlainText() #dpg
        sigmag = self.Inputsegmag.toPlainText() #sigma_g
        InitNumbConc = self.InputInitNumbConc.toPlainText()
        CoagFlag = self.InputCoag.toPlainText()
        CondFlag = self.InputCond.toPlainText()
        NuclFlag = self.InputNucl.toPlainText()
        ReacFlag = self.InputReac.toPlainText()
        TotT = self.InputFinalSimuTIme.toPlainText()
        ReacRate = self.InputReactRate.toPlainText()
                
        #OutputonScreen = Temp + '\n' + MW + '\n' + Den + '\n' + FD + '\n' + InitNumbConc + '\n' \
        #                 + DisNumb + '\n' + CoagFlag + '\n' + TotT + '\n' + Dt
        
        #self.OutputScreen.setText(OutputonScreen)
        
        outfile = open('input1.txt', 'w')
        
        output = InitNumbConc + ' ' + dpg + ' ' + sigmag + '\n' \
                 + SatPressure + '\n' \
                 + MW + ' ' + Den + ' ' + SurfT + '\n' \
                 + CoagFlag + '\n' \
                 + CondFlag + '\n' \
                 + NuclFlag + '\n' \
                 + ReacFlag + '\n' \
                 + ReacRate 
        
        outfile.write(output)
        
        outfile.close()
        
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Moment_Dialog = QtWidgets.QDialog()
    ui = Ui_Moment_Dialog()
    ui.setupUi(Moment_Dialog)
    Moment_Dialog.show()
    sys.exit(app.exec_())

