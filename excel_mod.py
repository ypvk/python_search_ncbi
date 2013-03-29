import xlwt

class Excel:

    def __init__(self, name=''):
        self.filename = name
        self.values = []
        self.heads = []
    
    def setHead(self, heads):
        self.heads = heads

    def addValues(self, values):
        self.values.append(values)
    
    def setValues(self, values):
        self.values = values

    def clearValues(self):
        self.filename = ''
        self.values = []
        self.heads = []

    def setFilename(self, filename):
        self.filename = filename

    def save(self):
        wbk = xlwt.Workbook()
        sheet = wbk.add_sheet('sheet 1')
        
        # write head
        for i, value in enumerate(self.heads):
            sheet.write(0, i, value)
        # write values
        for i, values in enumerate(self.values):
            for j in range(len(self.heads)):
                value = ''
                if values.has_key(self.heads[j]):
                    value = values[self.heads[j]]
                sheet.write(i+1, j, '%s' % value)
        wbk.save(self.filename)

