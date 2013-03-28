from PySide import QtGui, QtCore
from PySide.QtGui import QApplication, QWidget, QPushButton, QLineEdit, QTextEdit, \
        QVBoxLayout, QHBoxLayout, QComboBox, QLabel
from PySide.QtCore import QThread

from Bio import Entrez

email = 'ypseu.2007@163.com'
MAX_THREAD = 10
RET_MAX = 20
all_output = []

class Dialog(QWidget):

    def __init__(self, parent=None):
        super(Dialog, self).__init__(parent)
        self.setupUI()
        self.setupConnection()
        self.threadPool = []
        #self.initSearchDb()

    def setupUI(self):

        self.pushButton = QPushButton(u"Search", self)
        self.lineEdit = QLineEdit(self)
        self.textEdit = QTextEdit(self)
        self.comboBox = QComboBox(self)
        self.label = QLabel(u"DB:", self)
        
        self.textEdit.setReadOnly(True)

        self.layout = QVBoxLayout()
        self.topLayout = QHBoxLayout()
        
        self.topLayout.addWidget(self.label)
        self.topLayout.addWidget(self.comboBox)
        self.topLayout.addWidget(self.lineEdit)
        self.topLayout.addWidget(self.pushButton)

        self.layout.addLayout(self.topLayout)
        self.layout.addWidget(self.textEdit)

        self.setLayout(self.layout)

        self.resize(600, 700)
        self.setWindowTitle(u"Search Data for NCBI")

    def setupConnection(self):
        self.pushButton.clicked.connect(self.onButtonClicked)

    def onButtonClicked(self):
        if not self.lineEdit.text() or not self.comboBox.count():
            QtGui.QMessageBox.information(self, u"Warning", u"Please Set the Search Field")
            return
        dbName = self.comboBox.currentText()
        fieldName = self.lineEdit.text()
        self.log("Start searching db: %s and field: %s" % (dbName, fieldName))
        handle = self.entrez.esearch(db=dbName, term=fieldName, rettype='count')
        record = self.entrez.read(handle)
        self.log("result %s" % record)
        #handle = self.entrez.efetch(db=dbName, id=record['IdList'], rettype='gb')
        self.finishedThreadNum = 0
        if int(record['Count']) / RET_MAX > MAX_THREAD:
            self.realThreadNum = MAX_THREAD
            each_count = int(record['Count'])/MAX_THREAD
            start = 0
            for i in range(MAX_THREAD - 1):
                thread = MyThread(start=start, count=each_count, dbName=dbName, fieldName=fieldName)
                thread.finished.connect(self.onThreadFinished)
                print thread
                print dir(thread)
                thread.start()
                print thread
                self.threadPool.append(thread)
                start = start + each_count
            thread = MyThread(start=start, count=(int(record['Count'])-start+1), dbName=dbName, fieldName=fieldName)
            thread.finished.connect(self.onThreadFinished)
            self.threadPool.append(thread)
            thread.start()
        else:
            if int(record['Count']) == RET_MAX:
                self.realThreadNum = 1
            else:
                self.realThreadNum = int(record['Count'])/RET_MAX + 1
            start = 0
            for i in range(self.realThreadNum):
                thread = MyThread(start=start, count=RET_MAX, dbName=dbName, fieldName=fieldName)
                thread.finished.connect(self.onThreadFinished)
                self.threadPool.append(thread)
                thread.start()
                start = start + RET_MAX
        self.log('reading data')

    def log(self, context):
        self.textEdit.append(context)

    def initSearchDb(self):
        self.entrez = Entrez 
        self.entrez.email = email
        self.log("Connect to NCBI")
        handle = self.entrez.einfo()
        record = self.entrez.read(handle)
        self.log("Get NCBI DataBases:")
        for name in record['DbList']:
            self.log('DBName:\t%s' % name)
        self.comboBox.addItems(record['DbList'])

    def onThreadFinished(self):
        self.finishedThreadNum = self.finishedThreadNum + 1
        self.log('finished %s ' % self.finished)
        if(self.finishedThreadNum == self.realThreadNum):
            print all_output

class MyThread(QThread):
    def __init__(self, start, count, dbName, fieldName):
        super(MyThread, self).__init__()
        #self.entrez = Entrez 
        #self.entrez.email = email
        #self.start = start
        #self.count = count
        #self.dbName = dbName
        #self.fieldName = fieldName
        #self.output = []

    def run(self):
        print "yuping"
        #times = None 
        #if self.count == RET_MAX:
            #times = 1
        #else:
            #times = self.count/RET_MAX + 1
        #n = 0
        #while n < times:
            #handle = self.entrez.esearch(db=self.dbName, term=self.fieldName, usehistory='y', retstart=(self.start + n*RET_MAX))
            #record = self.entrez.read(handle)
            #handle = self.entrez.efetch(db=dbName, id=record['IdList'], rettype='gb')
            #self.output.append(handle.read())
            #n = n+1
        #all_output.append(self.output)
        self.exec_() 
