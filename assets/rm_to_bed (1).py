#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as rpmsk:
    with open(sys.argv[2],'w') as bed:
        for line in rpmsk:
            spl=line.split()
            if spl[0][:1].isdigit():
                if len(spl) > 5:
                    orient=spl[8]
                    if orient == 'C':
                        bed.write(spl[4]+'\t'+spl[5]+'\t'+spl[6]+'\t'+spl[9]+'#'+spl[10]+'\t'+'.'+'\t'+'-'+'\n')
                    else:
                        bed.write(spl[4]+'\t'+spl[5]+'\t'+spl[6]+'\t'+spl[9]+'#'+spl[10]+'\t'+'.'+'\t'+'+'+'\n')
