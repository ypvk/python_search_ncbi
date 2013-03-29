from excel_mod import Excel

excel = Excel("test.xls")
heads = ['a', 'n']
values = [{'a': 1, 'n':2}, {'n': 3, 'a': 10}]
excel.setHead(heads)
excel.setValues(values)
excel.save()
