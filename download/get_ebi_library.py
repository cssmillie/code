import re, urllib2

fail = open('failures.txt', 'w')

for sid in open('ids.lst'):
    sid = sid.rstrip()
    url = 'http://www.ebi.ac.uk/ena/data/view/%s&display=xml' %(sid)
    try:
        html = urllib2.urlopen(url).read()
        lib = re.search('qiita_ptid_\d+\:(.*?)\"', html).group(1)
        print sid, lib
    except:
        fail.write('%s\n' %(sid))
fail.close()
