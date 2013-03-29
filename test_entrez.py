import sys
from Bio import Entrez

if __name__ == "__main__":
    Entrez.email = "test@163.com"
    dbName = sys.argv[1]
    fieldName = sys.argv[2]
    handle = Entrez.esearch(db=dbName, term=fieldName, usehistory='y')
    record = Entrez.read(handle)
    count = int(record["Count"])
    WebEnv = record['WebEnv']
    QueryKey = record['QueryKey']
    print "%s_%s_%s" % (count, WebEnv, QueryKey)
