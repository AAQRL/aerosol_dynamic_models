# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\huangzhang\Desktop\Aerosol dynamics modeling\GUI\PyQt5 GUI\ADM\Dis_Dis-Sec_Moment_Modal\Modal_Dialog.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Modal_Dialog(object):
    
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
    
    def setupUi(self, Modal_Dialog):
        Modal_Dialog.setObjectName("Modal_Dialog")
        Modal_Dialog.resize(543, 466)
        self.COMSOLHere = QtWidgets.QLabel(Modal_Dialog)
        self.COMSOLHere.setGeometry(QtCore.QRect(130, 260, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.COMSOLHere.setFont(font)
        self.COMSOLHere.setObjectName("COMSOLHere")
        self.InputFinalSimuTIme = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputFinalSimuTIme.setGeometry(QtCore.QRect(280, 310, 104, 21))
        self.InputFinalSimuTIme.setObjectName("InputFinalSimuTIme")
        self.InputCond = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputCond.setGeometry(QtCore.QRect(280, 200, 104, 21))
        self.InputCond.setObjectName("InputCond")
        self.InputCOMSOL = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputCOMSOL.setGeometry(QtCore.QRect(280, 270, 104, 21))
        self.InputCOMSOL.setObjectName("InputCOMSOL")
        self.InputNucl = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputNucl.setGeometry(QtCore.QRect(280, 240, 104, 21))
        self.InputNucl.setObjectName("InputNucl")
        self.InputInitVaporPressure = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputInitVaporPressure.setGeometry(QtCore.QRect(280, 90, 104, 21))
        self.InputInitVaporPressure.setObjectName("InputInitVaporPressure")
        self.InputVaporGeneRate = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputVaporGeneRate.setGeometry(QtCore.QRect(280, 130, 104, 21))
        self.InputVaporGeneRate.setObjectName("InputVaporGeneRate")
        self.InitVaporPressureHere = QtWidgets.QLabel(Modal_Dialog)
        self.InitVaporPressureHere.setGeometry(QtCore.QRect(110, 80, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.InitVaporPressureHere.setFont(font)
        self.InitVaporPressureHere.setObjectName("InitVaporPressureHere")
        self.FinalSimuTimeHere = QtWidgets.QLabel(Modal_Dialog)
        self.FinalSimuTimeHere.setGeometry(QtCore.QRect(140, 300, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.FinalSimuTimeHere.setFont(font)
        self.FinalSimuTimeHere.setObjectName("FinalSimuTimeHere")
        self.CondHere = QtWidgets.QLabel(Modal_Dialog)
        self.CondHere.setGeometry(QtCore.QRect(100, 190, 181, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CondHere.setFont(font)
        self.CondHere.setObjectName("CondHere")
        self.NumbModesHere = QtWidgets.QLabel(Modal_Dialog)
        self.NumbModesHere.setGeometry(QtCore.QRect(150, 40, 111, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NumbModesHere.setFont(font)
        self.NumbModesHere.setObjectName("NumbModesHere")
        self.InputNumbModes = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputNumbModes.setGeometry(QtCore.QRect(280, 50, 104, 21))
        self.InputNumbModes.setObjectName("InputNumbModes")
        self.InputCoag = QtWidgets.QTextEdit(Modal_Dialog)
        self.InputCoag.setGeometry(QtCore.QRect(280, 160, 104, 21))
        self.InputCoag.setObjectName("InputCoag")
        self.CoagHere = QtWidgets.QLabel(Modal_Dialog)
        self.CoagHere.setGeometry(QtCore.QRect(110, 150, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.CoagHere.setFont(font)
        self.CoagHere.setObjectName("CoagHere")
        self.VaporGeneRateHere = QtWidgets.QLabel(Modal_Dialog)
        self.VaporGeneRateHere.setGeometry(QtCore.QRect(70, 120, 201, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.VaporGeneRateHere.setFont(font)
        self.VaporGeneRateHere.setObjectName("VaporGeneRateHere")
        self.NuclHere = QtWidgets.QLabel(Modal_Dialog)
        self.NuclHere.setGeometry(QtCore.QRect(120, 230, 161, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.NuclHere.setFont(font)
        self.NuclHere.setObjectName("NuclHere")
        self.RunModal = QtWidgets.QPushButton(Modal_Dialog)
        self.RunModal.setGeometry(QtCore.QRect(210, 380, 75, 23))
        self.RunModal.setObjectName("RunModal")

        self.retranslateUi(Modal_Dialog)
        QtCore.QMetaObject.connectSlotsByName(Modal_Dialog)
        
        #HZ 2019
        self.RunModal.clicked.connect(self.get_Modal_inputfile)
        #HZ 2019

    def retranslateUi(self, Modal_Dialog):
        _translate = QtCore.QCoreApplication.translate
        Modal_Dialog.setWindowTitle(_translate("Modal_Dialog", "Dialog"))
        self.COMSOLHere.setText(_translate("Modal_Dialog", "COMSOL Flag (0 or 1)"))
        self.InputFinalSimuTIme.setHtml(_translate("Modal_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCond.setHtml(_translate("Modal_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputCOMSOL.setHtml(_translate("Modal_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InputNucl.setHtml(_translate("Modal_Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.InitVaporPressureHere.setText(_translate("Modal_Dialog", "Initial Vapor Pressure (Pa)"))
        self.FinalSimuTimeHere.setText(_translate("Modal_Dialog", "Final Simu. Time (s)"))
        self.CondHere.setText(_translate("Modal_Dialog", "Condensation Flag (0 or 1)"))
        self.NumbModesHere.setText(_translate("Modal_Dialog", "Number of Modes"))
        self.CoagHere.setText(_translate("Modal_Dialog", "Coagulation  Flag (0 or 1)"))
        self.VaporGeneRateHere.setText(_translate("Modal_Dialog", "Vapor Generation Rate (#/m^3/s)"))
        self.NuclHere.setText(_translate("Modal_Dialog", "Nucleation Flag (0 or 1)"))
        self.RunModal.setText(_translate("Modal_Dialog", "Run"))
        
    def get_Modal_inputfile(self):
        
        global Temp, Pressure, SatPressure, Den, MW, SurfT
        print("Temp  = " + Temp)
        print("Pressure  = " + Pressure)
        print("SatPressure  = " + SatPressure)
        print("MW  = " + MW)
        print("Den  = " + Den)
        print("SurfT  = " + SurfT)
        
        NumbModes = self.InputNumbModes.toPlainText()
        InitVaporPressure = self.InputInitVaporPressure.toPlainText()
        VaporGeneRate = self.InputVaporGeneRate.toPlainText()
        CoagFlag = self.InputCoag.toPlainText()
        CondFlag = self.InputCond.toPlainText()
        NuclFlag = self.InputNucl.toPlainText()
        ComsolFlag = self.InputCOMSOL.toPlainText()
        TotT = self.InputFinalSimuTIme.toPlainText()
    
        outfile = open('input_file.txt', 'w')
        output = str(1) + '\n' \
                 + InitVaporPressure + '\n' \
                 + SatPressure + '\n' \
                 + VaporGeneRate + '\n' \
                 + '\n' \
                 + str(2) + '\n' \
                 + MW + '\n' \
                 + Den + '\n' \
                 + SurfT + '\n' \
                 + '\n' \
                 + str(3) + '\n' \
                 + str(29) + '\n' \
                 + str(1.255) + '\n' \
                 + str(1.9533e-5) + '\n' \
                 + str(110.4) + '\n' \
                 + str(67.4e-9) + '\n' \
                 + '\n' \
                 + str(4) + '\n' \
                 + Pressure+ '\n' \
                 + '\n' \
                 + str(5) + '\n' \
                 + NumbModes + '\n' \
                 + '\n' \
                 + str(6) + '\n' \
                 + CoagFlag + '\n' \
                 + CondFlag + '\n' \
                 + NuclFlag + '\n' \
                 + str(1) 
                 
        outfile.write(output)
        outfile.close()
        
        outfile = open('streamline_hist.txt', 'w')
        
        output  = str(0) + ' ' + str(1) + ' ' + Temp + '\n' \
                  + str(15) + ' ' + str(1) + ' ' + Temp
        
        outfile.write(output)
        outfile.close()
        
        NumbModes = int(NumbModes)
        #print('NumbModels', NumbModels)
        
        outfile1 = open('input_N.txt', 'w')
        outfile2 = open('input_v.txt', 'w')
        
        for i in range(NumbModes):
            output1 = str(1e-40) + '\n'
            outfile1.write(output1)
            output2 = str(0) + '\n'
            outfile2.write(output2)
        
        outfile1.close()
        outfile2.close()        
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Modal_Dialog = QtWidgets.QDialog()
    ui = Ui_Modal_Dialog()
    ui.setupUi(Modal_Dialog)
    Modal_Dialog.show()
    sys.exit(app.exec_())

