# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI_1D.ui'
#
# Created: Fri Mar 21 14:42:31 2008
#      by: PyQt4 UI code generator 4.3.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(QtCore.QSize(QtCore.QRect(0,0,768,605).size()).expandedTo(Dialog.minimumSizeHint()))

        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(520,540,171,35))

        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.buttonBox.setFont(font)
        self.buttonBox.setCursor(QtCore.Qt.ArrowCursor)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.NoButton|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")

        self.InFileGroupBox = QtGui.QGroupBox(Dialog)
        self.InFileGroupBox.setGeometry(QtCore.QRect(10,40,461,121))
        self.InFileGroupBox.setObjectName("InFileGroupBox")

        self.widget = QtGui.QWidget(self.InFileGroupBox)
        self.widget.setGeometry(QtCore.QRect(10,30,441,78))
        self.widget.setObjectName("widget")

        self.gridlayout = QtGui.QGridLayout(self.widget)
        self.gridlayout.setObjectName("gridlayout")

        self.lineEdit = QtGui.QLineEdit(self.widget)
        self.lineEdit.setObjectName("lineEdit")
        self.gridlayout.addWidget(self.lineEdit,0,0,1,1)

        self.GetParamFileButton = QtGui.QPushButton(self.widget)
        self.GetParamFileButton.setObjectName("GetParamFileButton")
        self.gridlayout.addWidget(self.GetParamFileButton,0,1,1,1)

        self.lineEdit_2 = QtGui.QLineEdit(self.widget)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.gridlayout.addWidget(self.lineEdit_2,1,0,1,1)

        self.GetSignalFileButton = QtGui.QPushButton(self.widget)
        self.GetSignalFileButton.setObjectName("GetSignalFileButton")
        self.gridlayout.addWidget(self.GetSignalFileButton,1,1,1,1)

        self.OutFileGroupBox = QtGui.QGroupBox(Dialog)
        self.OutFileGroupBox.setGeometry(QtCore.QRect(10,170,481,331))
        self.OutFileGroupBox.setObjectName("OutFileGroupBox")

        self.lineEdit_3 = QtGui.QLineEdit(self.OutFileGroupBox)
        self.lineEdit_3.setGeometry(QtCore.QRect(21,27,306,29))
        self.lineEdit_3.setObjectName("lineEdit_3")

        self.FDMLineListButton = QtGui.QPushButton(self.OutFileGroupBox)
        self.FDMLineListButton.setGeometry(QtCore.QRect(333,24,127,35))
        self.FDMLineListButton.setObjectName("FDMLineListButton")

        self.lineEdit_4 = QtGui.QLineEdit(self.OutFileGroupBox)
        self.lineEdit_4.setGeometry(QtCore.QRect(21,159,306,29))
        self.lineEdit_4.setObjectName("lineEdit_4")

        self.AbsSpButton = QtGui.QPushButton(self.OutFileGroupBox)
        self.AbsSpButton.setGeometry(QtCore.QRect(333,156,127,35))
        self.AbsSpButton.setObjectName("AbsSpButton")

        self.lineEdit_5 = QtGui.QLineEdit(self.OutFileGroupBox)
        self.lineEdit_5.setGeometry(QtCore.QRect(21,71,306,29))
        self.lineEdit_5.setObjectName("lineEdit_5")

        self.ReSpButton = QtGui.QPushButton(self.OutFileGroupBox)
        self.ReSpButton.setGeometry(QtCore.QRect(333,68,127,35))
        self.ReSpButton.setObjectName("ReSpButton")

        self.lineEdit_6 = QtGui.QLineEdit(self.OutFileGroupBox)
        self.lineEdit_6.setGeometry(QtCore.QRect(21,115,306,29))
        self.lineEdit_6.setObjectName("lineEdit_6")

        self.ImSpButton = QtGui.QPushButton(self.OutFileGroupBox)
        self.ImSpButton.setGeometry(QtCore.QRect(333,112,127,35))
        self.ImSpButton.setObjectName("ImSpButton")

        self.groupBox = QtGui.QGroupBox(self.OutFileGroupBox)
        self.groupBox.setGeometry(QtCore.QRect(10,190,461,131))
        self.groupBox.setObjectName("groupBox")

        self.radioButton = QtGui.QRadioButton(self.groupBox)
        self.radioButton.setGeometry(QtCore.QRect(10,100,127,28))
        self.radioButton.setObjectName("radioButton")

        self.label_12 = QtGui.QLabel(self.groupBox)
        self.label_12.setGeometry(QtCore.QRect(70,70,78,24))
        self.label_12.setObjectName("label_12")

        self.lineEdit_17 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_17.setGeometry(QtCore.QRect(10,70,51,29))
        self.lineEdit_17.setObjectName("lineEdit_17")

        self.lineEdit_10 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_10.setGeometry(QtCore.QRect(0,30,306,29))
        self.lineEdit_10.setObjectName("lineEdit_10")

        self.ParButton = QtGui.QPushButton(self.groupBox)
        self.ParButton.setGeometry(QtCore.QRect(320,30,127,35))
        self.ParButton.setObjectName("ParButton")

        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(330,70,101,24))
        self.label_6.setObjectName("label_6")

        self.lineEdit_11 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_11.setGeometry(QtCore.QRect(230,70,81,29))
        self.lineEdit_11.setObjectName("lineEdit_11")

        self.thetaLineEdit = QtGui.QLineEdit(Dialog)
        self.thetaLineEdit.setGeometry(QtCore.QRect(490,130,61,29))
        self.thetaLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.thetaLineEdit.setObjectName("thetaLineEdit")

        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(560,130,51,24))
        self.label.setObjectName("label")

        self.comboBox = QtGui.QComboBox(Dialog)
        self.comboBox.setGeometry(QtCore.QRect(490,50,71,28))
        self.comboBox.setObjectName("comboBox")

        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(590,60,78,24))
        self.label_2.setObjectName("label_2")

        self.lineEdit_7 = QtGui.QLineEdit(Dialog)
        self.lineEdit_7.setGeometry(QtCore.QRect(490,90,81,29))
        self.lineEdit_7.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_7.setObjectName("lineEdit_7")

        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(590,90,41,24))
        self.label_3.setObjectName("label_3")

        self.lineEdit_12 = QtGui.QLineEdit(Dialog)
        self.lineEdit_12.setGeometry(QtCore.QRect(630,130,71,29))
        self.lineEdit_12.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_12.setObjectName("lineEdit_12")

        self.label_7 = QtGui.QLabel(Dialog)
        self.label_7.setGeometry(QtCore.QRect(710,130,41,24))
        self.label_7.setObjectName("label_7")

        self.lineEdit_16 = QtGui.QLineEdit(Dialog)
        self.lineEdit_16.setGeometry(QtCore.QRect(540,480,113,29))
        self.lineEdit_16.setObjectName("lineEdit_16")

        self.label_11 = QtGui.QLabel(Dialog)
        self.label_11.setGeometry(QtCore.QRect(660,480,61,24))
        self.label_11.setObjectName("label_11")

        self.groupBox_2 = QtGui.QGroupBox(Dialog)
        self.groupBox_2.setGeometry(QtCore.QRect(10,510,211,80))
        self.groupBox_2.setObjectName("groupBox_2")

        self.lineEdit_18 = QtGui.QLineEdit(self.groupBox_2)
        self.lineEdit_18.setGeometry(QtCore.QRect(30,30,113,29))
        self.lineEdit_18.setObjectName("lineEdit_18")

        self.label_13 = QtGui.QLabel(self.groupBox_2)
        self.label_13.setGeometry(QtCore.QRect(160,30,78,24))
        self.label_13.setObjectName("label_13")

        self.groupBox_3 = QtGui.QGroupBox(Dialog)
        self.groupBox_3.setGeometry(QtCore.QRect(510,310,211,151))
        self.groupBox_3.setObjectName("groupBox_3")

        self.lineEdit_15 = QtGui.QLineEdit(self.groupBox_3)
        self.lineEdit_15.setGeometry(QtCore.QRect(10,110,71,29))
        self.lineEdit_15.setObjectName("lineEdit_15")

        self.lineEdit_13 = QtGui.QLineEdit(self.groupBox_3)
        self.lineEdit_13.setGeometry(QtCore.QRect(10,30,71,29))
        self.lineEdit_13.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_13.setObjectName("lineEdit_13")

        self.lineEdit_14 = QtGui.QLineEdit(self.groupBox_3)
        self.lineEdit_14.setGeometry(QtCore.QRect(10,70,71,29))
        self.lineEdit_14.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_14.setObjectName("lineEdit_14")

        self.label_10 = QtGui.QLabel(self.groupBox_3)
        self.label_10.setGeometry(QtCore.QRect(100,110,41,24))
        self.label_10.setObjectName("label_10")

        self.label_8 = QtGui.QLabel(self.groupBox_3)
        self.label_8.setGeometry(QtCore.QRect(100,30,51,24))
        self.label_8.setObjectName("label_8")

        self.label_9 = QtGui.QLabel(self.groupBox_3)
        self.label_9.setGeometry(QtCore.QRect(100,70,78,24))
        self.label_9.setObjectName("label_9")

        self.groupBox_4 = QtGui.QGroupBox(Dialog)
        self.groupBox_4.setGeometry(QtCore.QRect(510,180,211,111))
        self.groupBox_4.setObjectName("groupBox_4")

        self.lineEdit_9 = QtGui.QLineEdit(self.groupBox_4)
        self.lineEdit_9.setGeometry(QtCore.QRect(10,70,71,29))
        self.lineEdit_9.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_9.setObjectName("lineEdit_9")

        self.label_5 = QtGui.QLabel(self.groupBox_4)
        self.label_5.setGeometry(QtCore.QRect(90,70,78,24))
        self.label_5.setObjectName("label_5")

        self.lineEdit_8 = QtGui.QLineEdit(self.groupBox_4)
        self.lineEdit_8.setGeometry(QtCore.QRect(10,30,71,29))
        self.lineEdit_8.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_8.setObjectName("lineEdit_8")

        self.label_4 = QtGui.QLabel(self.groupBox_4)
        self.label_4.setGeometry(QtCore.QRect(90,30,78,24))
        self.label_4.setObjectName("label_4")

        self.groupBox_5 = QtGui.QGroupBox(Dialog)
        self.groupBox_5.setGeometry(QtCore.QRect(230,510,261,81))
        self.groupBox_5.setObjectName("groupBox_5")

        self.HelpButton = QtGui.QPushButton(self.groupBox_5)
        self.HelpButton.setGeometry(QtCore.QRect(10,30,111,35))
        self.HelpButton.setObjectName("HelpButton")

        self.programHelp = QtGui.QPushButton(self.groupBox_5)
        self.programHelp.setGeometry(QtCore.QRect(130,30,111,35))
        self.programHelp.setObjectName("programHelp")

        self.label_14 = QtGui.QLabel(Dialog)
        self.label_14.setGeometry(QtCore.QRect(310,0,131,40))

        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(20)
        font.setWeight(75)
        font.setBold(True)
        self.label_14.setFont(font)
        self.label_14.setTextFormat(QtCore.Qt.PlainText)
        self.label_14.setObjectName("label_14")

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox,QtCore.SIGNAL("accepted()"),Dialog.accept)
        QtCore.QObject.connect(self.buttonBox,QtCore.SIGNAL("rejected()"),Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.InFileGroupBox.setTitle(QtGui.QApplication.translate("Dialog", "Input Files", None, QtGui.QApplication.UnicodeUTF8))
        self.GetParamFileButton.setText(QtGui.QApplication.translate("Dialog", "Parameter File", None, QtGui.QApplication.UnicodeUTF8))
        self.GetSignalFileButton.setText(QtGui.QApplication.translate("Dialog", "Signal file", None, QtGui.QApplication.UnicodeUTF8))
        self.OutFileGroupBox.setTitle(QtGui.QApplication.translate("Dialog", "Output files", None, QtGui.QApplication.UnicodeUTF8))
        self.FDMLineListButton.setText(QtGui.QApplication.translate("Dialog", "FDM Linelist", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_4.setToolTip(QtGui.QApplication.translate("Dialog", "file where Abs[F(w)] will be written", None, QtGui.QApplication.UnicodeUTF8))
        self.AbsSpButton.setText(QtGui.QApplication.translate("Dialog", "Abs Spectrum", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_5.setToolTip(QtGui.QApplication.translate("Dialog", "file where Re[F(w)] will be written", None, QtGui.QApplication.UnicodeUTF8))
        self.ReSpButton.setText(QtGui.QApplication.translate("Dialog", "Re Spectrum", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_6.setToolTip(QtGui.QApplication.translate("Dialog", "file where Im[F(w)] will be written", None, QtGui.QApplication.UnicodeUTF8))
        self.ImSpButton.setText(QtGui.QApplication.translate("Dialog", "Im Spectrum", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("Dialog", "FDM Only", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton.setToolTip(QtGui.QApplication.translate("Dialog", "if checked, F(w) is computed with Im d_k (FDM only)", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton.setText(QtGui.QApplication.translate("Dialog", "CheatMore", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("Dialog", "Cheat", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_17.setToolTip(QtGui.QApplication.translate("Dialog", "multiply all widths by cheat (FDM only)", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_17.setText(QtGui.QApplication.translate("Dialog", "1", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_10.setToolTip(QtGui.QApplication.translate("Dialog", "file where the linelist will be written (FDM only)", None, QtGui.QApplication.UnicodeUTF8))
        self.ParButton.setText(QtGui.QApplication.translate("Dialog", "par ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("Dialog", "threshhold", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_11.setText(QtGui.QApplication.translate("Dialog", "0.0001", None, QtGui.QApplication.UnicodeUTF8))
        self.thetaLineEdit.setToolTip(QtGui.QApplication.translate("Dialog", "zero order (overall) phase correction", None, QtGui.QApplication.UnicodeUTF8))
        self.thetaLineEdit.setText(QtGui.QApplication.translate("Dialog", "1.5", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Dialog", "theta", None, QtGui.QApplication.UnicodeUTF8))
        self.comboBox.setToolTip(QtGui.QApplication.translate("Dialog", "method used to compute the spectra", None, QtGui.QApplication.UnicodeUTF8))
        self.comboBox.addItem(QtGui.QApplication.translate("Dialog", "FDM", None, QtGui.QApplication.UnicodeUTF8))
        self.comboBox.addItem(QtGui.QApplication.translate("Dialog", "RRT", None, QtGui.QApplication.UnicodeUTF8))
        self.comboBox.addItem(QtGui.QApplication.translate("Dialog", "DFT", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Dialog", "Method", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_7.setText(QtGui.QApplication.translate("Dialog", "40960", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Dialog", "Nsig", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_12.setToolTip(QtGui.QApplication.translate("Dialog", " basis density", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_12.setText(QtGui.QApplication.translate("Dialog", "1.1", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("Dialog", "rho", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_16.setToolTip(QtGui.QApplication.translate("Dialog", "smoothing parameter. Defines the smallest possible linewidth.", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_16.setText(QtGui.QApplication.translate("Dialog", "5E-02", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("Dialog", "Gamm", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setTitle(QtGui.QApplication.translate("Dialog", "RRT Only", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_18.setToolTip(QtGui.QApplication.translate("Dialog", "regularization parameter (RRT only)", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_18.setText(QtGui.QApplication.translate("Dialog", "1e-08 ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("Dialog", "ros", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_3.setTitle(QtGui.QApplication.translate("Dialog", "Fourier Basis Options", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_15.setToolTip(QtGui.QApplication.translate("Dialog", "number of points used in plotting the spectra", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_15.setText(QtGui.QApplication.translate("Dialog", "20000", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_13.setToolTip(QtGui.QApplication.translate("Dialog", "number of narrow band Fourier basis functions used per window", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_13.setText(QtGui.QApplication.translate("Dialog", "1050", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_14.setText(QtGui.QApplication.translate("Dialog", "20", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("Dialog", "Nsp", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setToolTip(QtGui.QApplication.translate("Dialog", "number of narrow band Fourier basis functions used per window", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("Dialog", "Nb0", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("Dialog", "Nbc", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_4.setTitle(QtGui.QApplication.translate("Dialog", "Window Params", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_9.setToolTip(QtGui.QApplication.translate("Dialog", "specify the mimimum value of the spectra to be computed", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_9.setText(QtGui.QApplication.translate("Dialog", "-9000", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("Dialog", "wmino", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_8.setToolTip(QtGui.QApplication.translate("Dialog", "maximum value of the spectra to be computed.", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_8.setText(QtGui.QApplication.translate("Dialog", "4500", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("Dialog", "wmaxo", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_5.setTitle(QtGui.QApplication.translate("Dialog", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.HelpButton.setText(QtGui.QApplication.translate("Dialog", "on Options", None, QtGui.QApplication.UnicodeUTF8))
        self.programHelp.setText(QtGui.QApplication.translate("Dialog", "on Program", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("Dialog", "fd_rr1d", None, QtGui.QApplication.UnicodeUTF8))

