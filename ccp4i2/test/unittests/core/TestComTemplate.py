import os
import unittest

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.CCP4Utils import getCCP4I2Dir
from ccp4i2.core.CCP4ComTemplate import CComTemplate, CComTemplateIf, CComTemplateLine


def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testComTemplate)
    #suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testQObject))
    return suite


def testModule():
    suite = TESTSUITE()
    #print 'suite',suite
    unittest.TextTestRunner(verbosity=2).run(suite)


class testComTemplate(unittest.TestCase):

    def makeContainer(self):
        c = CContainer()
        c.loadContentsFromXml(os.path.join(getCCP4I2Dir(),'test','data','test_com_template_1.def.xml'))
        c.loadDataFromXml(os.path.join(getCCP4I2Dir(),'test','data','test_com_template_1.params.xml'))
        #print 'makeContainer',c.dataOrder()
        return c

    def test_1(self):
        template = '''
1 whatever $WHATEVER
$TEST_PAR test $TEST_VALUE
'''
        t = CComTemplate(template=template)
        #print t
        self.assertEqual(len(t),2,'Failed to load template')

    def test_2(self):
        template = '''
1 title $TITLE
IF $HOWEVER
  1 whatever $WHATEVER
ELSE
  $TEST_PAR test $TEST_VALUE
ENDIF
1 trailer $TRAILER
'''
        t = CComTemplate(template=template)
        #print t
        self.assertEqual(len(t),3,'Failed to load template')
        self.assertTrue(isinstance(t.contents[1],CComTemplateIf),'Failed to load template CComTemplateIf')

    def test_3(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v = t.getValue('$WHATEVER',container=c)
        self.assertEqual(v,'12','Failed to convert $WHATEVER to data value')
        v = t.getValue('$TITLE.text',container=c)
        self.assertEqual(v,'Try this title','Failed to convert $TITLE.text to data value')
        try: 
            v = t.getValue('$DONTEXIST',container=c)
        except CException as e:
            errcode = e[0]['code']
        self.assertEqual(errcode,1,'Failed to convert $DONTEXIST to correct error code')


    def test_4(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v,err = t.substituteValues('1 whatever $WHATEVER',container=c)
        self.assertEqual(v,'1 whatever 12','Failed to substituteValues in WHATEVER line')    
        v,err = t.substituteValues('$TEST_PAR test $TEST_VALUE',container=c)
        self.assertEqual(v,'True test 24.0','Failed to substituteValues in TEST_VALUE line')

    def test_5(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v,err = t.evaluateCodeToBoolean('$WHATEVER.isSet()',container=c)
        self.assertTrue(v,'Failed to evaluate $WHATEVER.isSet()')
        v,err = t.evaluateCodeToBoolean('10 < $WHATEVER < 15',container=c)
        self.assertTrue(v,'Failed to evaluate 10 < $WHATEVER < 15')

    def test_6(self):
        c = self.makeContainer()
        template = '''1 whatever $WHATEVER\n$TEST_PAR test $TEST_VALUE\n{$WHATEVER >20} greater $WHATEVER\n{$WHATEVER <20} less than $WHATEVER\n'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_6',text,err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\ntest 24.0\nless than 12\n''','Failed to create correct script')

    def test_7(self):
        c = self.makeContainer()
        template = '''1 whatever $WHATEVER\n- $TEST_PAR test $TEST_VALUE\n-- {$WHATEVER >20} greater $WHATEVER\n-- {$WHATEVER <20} less than $WHATEVER\n'''
        t = CComTemplate(template=template)
        #print 'test_7'
        #print t
        text,err = t.makeComScript(container=c)
        #print 'test_7'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\n  - test 24.0\n  - less than 12\n''','Failed to create correct script')

    def test_8(self):
        c = self.makeContainer()
        template = '''
1 whatever $WHATEVER
IF not $HOWEVER
  $TEST_PAR test $TEST_VALUE
  {$WHATEVER >20} greater $WHATEVER
ELSE
  {$WHATEVER <20} less than $WHATEVER\n
ENDIF
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_8'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\nless than 12\nfoo\n''','Failed to create correct script')

    def test_9(self):
        c = self.makeContainer()
        template = '''
1 whatever $WHATEVER
IF $HOWEVER
  IF $TEST_PAR
    1 test $TEST_VALUE
  ELSE
    1 greater $WHATEVER
  ENDIF
ELSE
  {$WHATEVER <20} less than $WHATEVER\n
ENDIF
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_9'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\ntest 24.0\nfoo\n''','Failed to create correct script')


    def test_10(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
LOOP N 1 $WHATEVER
  1 RUN $N
ENDLOOP
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text, err = t.makeComScript(container=c)
        #print 'test_10'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_11(self):
        l = CComTemplateLine()
        t,err = l.evaluateCodeFragments('test { 3 * 4 } this',None)
        self.assertEqual(t,'test 12 this','Error in CComTemplateLine.evaluateCodeFragments')
        self.assertEqual(len(err),0,'Unexpected errors in test 11')
 
    def test_12(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
LOOP N 1 $WHATEVER
  1 RUN $N
  - 1 TEST { $N * $TEST_VALUE }
ENDLOOP
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_12'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
 
    def test_13(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $TRAILER
CASEMATCH bar
  1 matched bar
CASEMATCH foo
  1 matched foo
ENDCASE
1 WHATEVER $WHATEVER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_13'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched foo\nWHATEVER 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_14(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $WHATEVER
CASEMATCH 12
  1 matched 12
CASEMATCH  24
  1 matched 24
ENDCASE
1 WHATEVER $WHATEVER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_14'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched 12\nWHATEVER 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_15(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $TRAILER
CASEMATCH bar
  1 matched bar
CASEMATCH foo
  1 matched foo
  AT $CCP4I2_TOP/test/data/test_com_template_1.com
ENDCASE
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_15'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched foo\nwhatever 12\nless than 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_16(self):
        c = self.makeContainer()
        template = '''
LABELLINE F $F_SIGF.F SIGF $F_SIGF.SIGF
LABELLINE PHI $PHI_W.PHI W $PHI_W.W
'''
        t = CComTemplate(template=template,diagnostic=False)
        text,err = t.makeComScript(container=c)
        #print 'test_16'
        #print text,text.split('\n')
        #print err
        self.assertEqual(len(err),1,'Wrong number of errors in LABELLINE test')
        self.assertEqual(text,'LABIN F=ffoo SIGF=ffoosig\nLABIN PHI=foophi\n','LABELLINE output incorrect')
