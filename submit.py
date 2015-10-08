#from sys import argv

#script, filename = argv

samplelist = ['WW50','Top50','TTJets50','HWW50','DY50','DY25','VV50','WJets50','Data201550']
channelist = ['OF','EE','MuMu','OF']
IDlist     = ['MediumIDTighterIP']
#['MediumID','MediumIDTighterIP','TightID','TightIDTighterIP']


for sample in samplelist:
    
    for channel in channelist:

        for ID in IDlist:

            filenames = "launch"+sample+channel+ID+".sh"
            import time
            time.sleep(1)
            print 'qsub %s' % filenames
            import os
            stream = os.system('qsub %s' % filenames)
            #stream = os.system('qstat -u trevisanin')
            
            
