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
    print "%s#%s#%s" % (count, WebEnv, QueryKey)
    RET_MAX_SUMMARY = 10000
    times = count / RET_MAX_SUMMARY + 1
    for i in range(times):
        handle = Entrez.esummary(db=dbName, retstart=i*RET_MAX_SUMMARY,
                retmax=RET_MAX_SUMMARY,
                WebEnv=WebEnv,
                query_key=QueryKey
                )
        result = Entrez.read(handle)
        print len(result)
        print result[0]
