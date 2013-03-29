from PySide import QtGui, QtCore
from PySide.QtGui import QApplication, QWidget, QPushButton, QLineEdit, QTextEdit, \
        QVBoxLayout, QHBoxLayout, QComboBox, QLabel, QProgressBar
from PySide.QtCore import QThread

from Bio import Entrez
from excel_mod import Excel

email = 'mytest@163.com'
MAX_THREAD = 10
RET_MAX = 20
all_output = []
RET_MAX_SUMMARY = 10000
MAX_COUNT = 100000
FINISHED_COUNT = 0

class Dialog(QWidget):

    def __init__(self, parent=None):
        super(Dialog, self).__init__(parent)
        self.setupUI()
        self.setupConnection()
        self.threadPool = []
        self.totalUids = 0
        self.finishedThreadNum = 0
        self.realThreadNum = 0
        self.excel = Excel()
        #self.initSearchDb()

    def setupUI(self):

        self.pushButton = QPushButton(u"Search", self)
        #self.testButton = QPushButton(u"Test", self)
        self.lineEdit = QLineEdit(self)
        self.textEdit = QTextEdit(self)
        self.comboBox = QComboBox(self)
        self.label = QLabel(u"DB:", self)
        self.progressBar = QProgressBar(self)
        self.progressBar.setRange(0,1)
        
        self.textEdit.setReadOnly(True)

        self.layout = QVBoxLayout()
        self.topLayout = QHBoxLayout()
        
        self.topLayout.addWidget(self.label)
        self.topLayout.addWidget(self.comboBox)
        self.topLayout.addWidget(self.lineEdit)
        self.topLayout.addWidget(self.pushButton)
        #self.topLayout.addWidget(self.testButton)

        #self.testButton.clicked.connect(self.onTestButtonClicked)

        self.layout.addLayout(self.topLayout)
        self.layout.addWidget(self.textEdit)
        self.layout.addWidget(self.progressBar)

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
        # add use history to add all uids to the history server
        handle = self.entrez.esearch(db=dbName, term=fieldName, usehistory='y')
        record = self.entrez.read(handle)
        self.log("All result count %s" % record['Count'])
        self.totalUids = record['Count']
        #handle = self.entrez.efetch(db=dbName, id=record['IdList'], rettype='gb')
        self.finishedThreadNum = 0
        WebEnv = record['WebEnv']
        QueryKey = record['QueryKey']
        FINISHED_COUNT = 0
        self.progressBar.setValue(0)
        self.progressBar.setMaximum(int(self.totalUids))
        if int(record['Count']) / RET_MAX_SUMMARY > MAX_THREAD:
            self.realThreadNum = MAX_THREAD
            each_count = int(record['Count'])/MAX_THREAD
            startIndex = 0
            for i in range(MAX_THREAD - 1):
                thread = MyThread(startIndex, each_count, dbName, fieldName, WebEnv, QueryKey)
                thread.finished.connect(self.onThreadFinished)
                thread.finishedCountChanged.connect(self.onFinishedCountChange)
                thread.start()
                self.threadPool.append(thread)
                startIndex = startIndex + each_count
            thread = MyThread(startIndex, (int(record['Count'])-startIndex+1), dbName, fieldName, WebEnv, QueryKey)
            thread.finished.connect(self.onThreadFinished)
            thread.finishedCountChanged.connect(self.onFinishedCountChange)
            self.threadPool.append(thread)
            thread.start()
        else:
            if int(record['Count']) == RET_MAX_SUMMARY:
                self.realThreadNum = 1
            else:
                self.realThreadNum = int(record['Count'])/RET_MAX_SUMMARY + 1
            startIndex = 0
            for i in range(self.realThreadNum):
                thread = MyThread(startIndex, RET_MAX_SUMMARY, dbName, fieldName, WebEnv, QueryKey)
                thread.finished.connect(self.onThreadFinished)
                thread.finishedCountChanged.connect(self.onFinishedCountChange)
                self.threadPool.append(thread)
                thread.start()
                startIndex = startIndex + RET_MAX_SUMMARY
        self.log("%s" % self.realThreadNum)
        self.log('reading data')
        filename = '%s_%s_output.xls' % (dbName, fieldName)
        self.excel.setFilename(filename)

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
        self.log('finished thread %s ' % self.finishedThreadNum)
        if(self.finishedThreadNum == self.realThreadNum):
            heads = all_output[0][0].keys()
            self.excel.setHead(heads)
            for values in all_output:
                for value in values:
                    self.excel.addValues(value)
            self.excel.save()
            self.progressBar.setValue(int(self.totalUids))
            # clear all thread
            self.threadPool = []
    def onFinishedCountChange(self, num):
        self.progressBar.setValue(num)

    def onTestButtonClicked(self):
        self.finishedThreadNum = 0
        self.realThreadNum = 1
        self.thread = MyThread(0, 10,"abd","123")
        self.thread.finished.connect(self.onThreadFinished)
        self.thread.startIndex()

class MyThread(QThread):
    # define signal
    finishedCountChanged = QtCore.Signal(int)

    def __init__(self, startIndex, count, dbName, fieldName, WebEnv, query_key):
        super(MyThread, self).__init__()
        self.entrez = Entrez 
        self.entrez.email = email
        self.startIndex = startIndex
        self.count = count
        self.dbName = dbName
        self.fieldName = fieldName
        self.WebEnv = WebEnv
        self.query_key = query_key

    def run(self):
        times = None 
        if self.count == RET_MAX_SUMMARY:
            times = 1
        else:
            times = self.count/RET_MAX_SUMMARY + 1
        n = 0
        while n < times:
            handle = self.entrez.esummary(
                    db=self.dbName,
                    retstart=(self.startIndex+n*RET_MAX_SUMMARY),
                    retmax=RET_MAX_SUMMARY,
                    WebEnv=self.WebEnv,
                    query_key=self.query_key
                    )
            result = self.entrez.read(handle)
            if not isinstance(result, list):
                result = [result]
            all_output.append(result)
            FINISHED_COUNT =+ len(result)
            self.finishedCountChanged.emit(FINISHED_COUNT)
            n = n+1
