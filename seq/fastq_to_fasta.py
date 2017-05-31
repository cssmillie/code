import sys

#!/usr/bin/env python
import util
q = util.iter_fsq(sys.argv[1])
for record in q:
    print '>%s\n%s' %(record[0][1:], record[1])
