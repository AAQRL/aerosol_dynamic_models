# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\huangzhang\Desktop\Aerosol dynamics modeling\GUI\PyQt5 GUI\ADM\Dis_Dialog_v2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dis_Dialog(object):
    
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
        #print("MW  = " + MW)
        #print("Den  = " + Den)
                
    #HZ 2019   
    
    def setupUi(self, Dis_Dialog):
        Dis_Dialog.setObjectName("Dis_Dialog")
        Dis_Dialog.resize(543, 466)
        self.COMSOLHere = QtWidgets.QLabel(Dis_Dialog)
        self.COMSOLHere.setGeometry(QtCore.QRect(110, 250, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.COMSOLHere.setFont(font)
        self.COMSOLHere.setObjectName("COMSOLHere")
        self.InputFinalSimuTIme = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputFinalSimuTIme.setGeometry(QtCore.QRect(280, 300, 104, 21))
        self.InputFinalSimuTIme.setObjectName("InputFinalSimuTIme")
        self.InputCond = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputCond.setGeometry(QtCore.QRect(280, 190, 104, 21))
        self.InputCond.setObjectName("InputCond")
        self.InputCOMSOL = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputCOMSOL.setGeometry(QtCore.QRect(280, 260, 104, 21))
        self.InputCOMSOL.setObjectName("InputCOMSOL")
        self.InputNucl = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputNucl.setGeometry(QtCore.QRect(280, 230, 104, 21))
        self.InputNucl.setObjectName("InputNucl")
        self.InputInitNumbConc = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputInitNumbConc.setGeometry(QtCore.QRect(280, 110, 104, 21))
        self.InputInitNumbConc.setObjectName("InputInitNumbConc")
        self.InputDt = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputDt.setGeometry(QtCore.QRect(280, 330, 104, 21))
        self.InputDt.setObjectName("InputDt")
        self.IniNumbConcHere = QtWidgets.QLabel(Dis_Dialog)
        self.IniNumbConcHere.setGeometry(QtCore.QRect(110, 100, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.IniNumbConcHere.setFont(font)
        self.IniNumbConcHere.setObjectName("IniNumbConcHere")
        self.FinalSimuTimeHere = QtWidgets.QLabel(Dis_Dialog)
        self.FinalSimuTimeHere.setGeometry(QtCore.QRect(110, 290, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FinalSimuTimeHere.setFont(font)
        self.FinalSimuTimeHere.setObjectName("FinalSimuTimeHere")
        self.CondHere = QtWidgets.QLabel(Dis_Dialog)
        self.CondHere.setGeometry(QtCore.QRect(110, 180, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CondHere.setFont(font)
        self.CondHere.setObjectName("CondHere")
        self.DisNumbHere = QtWidgets.QLabel(Dis_Dialog)
        self.DisNumbHere.setGeometry(QtCore.QRect(110, 60, 111, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DisNumbHere.setFont(font)
        self.DisNumbHere.setObjectName("DisNumbHere")
        self.InputDisNumb = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputDisNumb.setGeometry(QtCore.QRect(280, 70, 104, 21))
        self.InputDisNumb.setObjectName("InputDisNumb")
        self.FirstDiamHere = QtWidgets.QLabel(Dis_Dialog)
        self.FirstDiamHere.setGeometry(QtCore.QRect(110, 30, 101, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FirstDiamHere.setFont(font)
        self.FirstDiamHere.setObjectName("FirstDiamHere")
        self.InputCoag = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputCoag.setGeometry(QtCore.QRect(280, 150, 104, 21))
        self.InputCoag.setObjectName("InputCoag")
        self.CoagHere = QtWidgets.QLabel(Dis_Dialog)
        self.CoagHere.setGeometry(QtCore.QRect(110, 140, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CoagHere.setFont(font)
        self.CoagHere.setObjectName("CoagHere")
        self.DtHere = QtWidgets.QLabel(Dis_Dialog)
        self.DtHere.setGeometry(QtCore.QRect(110, 320, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.DtHere.setFont(font)
        self.DtHere.setObjectName("DtHere")
        self.NuclHere = QtWidgets.QLabel(Dis_Dialog)
        self.NuclHere.setGeometry(QtCore.QRect(110, 220, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NuclHere.setFont(font)
        self.NuclHere.setObjectName("NuclHere")
        self.InputFirstDiam = QtWidgets.QTextEdit(Dis_Dialog)
        self.InputFirstDiam.setGeometry(QtCore.QRect(280, 40, 104, 21))
        self.InputFirstDiam.setObjectName("InputFirstDiam")
        self.RunDis = QtWidgets.QPushButton(Dis_Dialog)
        self.RunDis.setGeometry(QtCore.QRect(210, 380, 75, 23))
        self.RunDis.setObjectName("RunDis")

        self.retranslateUi(Dis_Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dis_Dialog)
        
        #HZ 2019
        self.RunDis.clicked.connect(self.get_dis_inputfile)
        #HZ 2019

    def retranslateUi(self, Dis_Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dis_Dialog.setWindowTitle(_translate("Dis_Dialog", "Discrete Dialog"))
        self.COMSOLHere.setText(_translate("Dis_Dialog", "COMSOL Flag (0 or 1)"))
        self.InputFinalSimuTIme.setHtml(_translate("Dis_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCond.setHtml(_translate("Dis_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCOMSOL.setHtml(_translate("Dis_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputNucl.setHtml(_translate("Dis_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.IniNumbConcHere.setText(_translate("Dis_Dialog", "Inital Numb Conc. (#/m^3)"))
        self.FinalSimuTimeHere.setText(_translate("Dis_Dialog", "Final Simu. Time (s)"))
        self.CondHere.setText(_translate("Dis_Dialog", "Condensation Flag (0 or 1)"))
        self.DisNumbHere.setText(_translate("Dis_Dialog", "Discrete Number"))
        self.FirstDiamHere.setText(_translate("Dis_Dialog", "First Diamer (m)"))
        self.CoagHere.setText(_translate("Dis_Dialog", "Coagulation  Flag (0 or 1)"))
        self.DtHere.setText(_translate("Dis_Dialog", "Dt (s)"))
        self.NuclHere.setText(_translate("Dis_Dialog", "Nucleation Flag (0 or 1)"))
        self.RunDis.setText(_translate("Dis_Dialog", "Run"))
        
    def get_dis_inputfile(self):
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
        
        FD = self.InputFirstDiam.toPlainText() #first discrete diameter
        InitNumbConc = self.InputInitNumbConc.toPlainText()
        DisNumb = self.InputDisNumb.toPlainText()
        CoagFlag = self.InputCoag.toPlainText()
        TotT = self.InputFinalSimuTIme.toPlainText()
        Dt = self.InputDt.toPlainText()
        
        #OutputonScreen = Temp + '\n' + MW + '\n' + Den + '\n' + FD + '\n' + InitNumbConc + '\n' \
        #                 + DisNumb + '\n' + CoagFlag + '\n' + TotT + '\n' + Dt
        
        #self.OutputScreen.setText(OutputonScreen)
        
        outfile = open('input.txt', 'w')
        
        output = CoagFlag + '\n' + DisNumb + '\n' + MW + ' ' + Den + ' ' + FD + \
                 '\n' + InitNumbConc + '\n' + TotT + '\n' + Dt + '\n' + str(1e6) 
        
        outfile.write(output)
        
        outfile.close()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dis_Dialog = QtWidgets.QDialog()
    ui = Ui_Dis_Dialog()
    ui.setupUi(Dis_Dialog)
    Dis_Dialog.show()
    sys.exit(app.exec_())

