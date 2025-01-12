#!/usr/bin/env python
##!/usr/local/Python-3.1.1/python

import sys


if __name__ == '__main__':
    file=sys.argv[1]
    ofile = sys.argv[2]
#    print('sys.path=', sys.path, '<br>')
    
    fp=open(file,"r")
    fo=open(ofile,"w")
    
    non_ascii = "àáâãäåæÀÁÂÃÄÅÆßç©¢ÇðÐèéêëÈÉÊË¡ìíîïÌÍÎÏñÑøòóôõöØÒÓÔÕÖÞþ®ùúûüµÙÚÛÜýÿÝ"
    ascii =     "aaaaaaaAAAAAAAbcccCdDeeeeEEEEiiiiiIIIInNooooooOOOOOOppRuuuuuUUUUyyY"
        
    for line in fp:
        for c in line:
            num=ord(c)
            if (num>127 and c in non_ascii):
                n=non_ascii.index('%s' %c)
                print(("Warning: %s (%d) is converted to %s (%d)" %(c,num,ascii[n],n)))
                c=ascii[n]
                
            elif num>127 and c not in non_ascii:
                print(("%s is non-ASCII (%d), repleced by '?'" %(c,num)))
                c='?'
                
            fo.write(c)

    fo.close() 
    fp.close()
            


    '''
#    s='^Cgð¨è¢µÜ·'
    n=0
    for x in s:
        print(n, x, len(x), ord(x) )
        n +=1
    '''
              
