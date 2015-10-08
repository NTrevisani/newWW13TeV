#from sys import argv

#script, filename = argv

samplelist = ['WW50','Top50','TTJets50','HWW50','DY50','DY25','VV50','WJets50','Data201550']
channelist = ['SF','EE','MuMu','OF']
IDlist     = ['MediumIDTighterIP']
#['MediumID','MediumIDTighterIP','TightID','TightIDTighterIP']

for sample in samplelist:
    
    for channel in channelist:

        for ID in IDlist:

            filenames = "launch"+sample+channel+ID+".sh"
            
            target = open(filenames, 'w')

            target.truncate()
            
            target.write("###################")
            target.write("\n")
            target.write("# FLAG definition #")
            target.write("\n")
            target.write("###################")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("# Request the Bourne Shell")
            target.write("\n")
            target.write("#$ -S /bin/bash")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("# Change the job name to 'hello_parallel'")
            target.write("\n")
            target.write("#$ -N hxx_9")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("# We are using the 'l.gaes' project")
            target.write("\n")
            target.write("#$ -P l.gaes")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("#if you want to receive notifications when your job starts/ends")
            target.write("\n")
            target.write("#$ -M nicolo.trevisani89@gmail.com")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("# Submit a parallel job of 16 slots")
            target.write("\n")
            target.write("#$ -pe mpi 8")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("#$ -cwd")
            target.write("\n")
            target.write("#$ -o output.out")
            target.write("\n")
            target.write("#$ -e output.err")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("# Make a reservation. Since the job is paralell and it might be that all the")
            target.write("\n")
            target.write("# resources are used, ensure that our resources are reserved.")
            target.write("\n")
            target.write("#$ -R y")
            target.write("\n")
            target.write("#################")
            target.write("\n")
            target.write("# Actual script #")
            target.write("\n")
            target.write("#################")
            target.write("\n")
            target.write("")
            target.write("\n")
            target.write("source /gpfs/csic_projects/cms/sw/ROOT/current/root/bin/thisroot.sh")
            target.write("\n")
            target.write("source /gpfs/csic_projects/cms/sw/PAF_releases/old/PAF_setup.sh")
            target.write("\n")
            target.write("\n")
            
            target.write("root -l -b -q \"RunMuonAnalyzer.C(\\\"%s\\\",\\\"%s\\\",\\\"OS\\\",\\\"Sequential\\\",0.04008,\\\"%s\\\")\"" % (sample, channel, ID))
            
            target.close()
            print 'ok'

import time
time.sleep(1)
            
for sample in samplelist:
    
    for channel in channelist:

        for ID in IDlist:

            import os
            stream = os.system('chmod +x launch*')#os.system('chmod +x %s' % filenames)
            #print 'chmod +x %s' % filenames 

#for sample in samplelist:
    
#    for channel in channelist:

#       for ID in IDlist:
            
            #           import time
            #           time.sleep(1)
            #           print 'qsub %s' % filenames
            #           stream = os.system('qsub %s' % filenames)
            #           stream = os.system('qstat -u trevisanin')
            
            
