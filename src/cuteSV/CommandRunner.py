from string import Template
import tempfile
import subprocess, signal, logging, os, stat, sys

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm
    
def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

def exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                            stderr=subprocess.STDOUT, close_fds=True,\
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout*60))  
    try:
        stdoutVal, stderrVal =  proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d" % (timeout)))
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return 214,None,None
    
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

class Command():
    def __init__(self, cmd, jobname, stdout, stderr):
        self.cmd = cmd
        self.jobname = jobname
        self.stdout = stdout
        self.stderr = stderr
    
    def asDict(self):
        return {"CMD":self.cmd, "JOBNAME":self.jobname, \
                "STDOUT":self.stdout, "STDERR":self.stderr}
    
class CommandRunner():
    """
    Uses a command template to run stuff. This is helpful for cluster commands
    and chunking several commands together
    """
    def __init__(self, template=None, njobs=0):
        """
        template: a string that will become the template for submitting to your cluster:
            #you can also go ahead and specify a string.Template
            default is to not submit to your cluster
            ${CMD} > ${STDOUT} 2> ${STDERR}
        njobs: (0)
            for clumping commands together and submitting them in a script
        """
        if template is None:
            template = "${CMD} > ${STDOUT} 2> ${STDERR}"
            self.runType = "Running"
        else:
            self.runType = "Submitting"
        self.template = Template(template)
        self.njobs = njobs
    
    def __call__(self, cmds, wDir = None, id = None):
        """
        Executes Commands - can either be a list or a single Command
        wDir is the working directory where chunk scripts will be written
        if id is None a random identifier will be applied when chunking
        """
        if wDir is None:
            wDir = "./"
        
        if type(cmds) != list:
            cmd = self.buildCommand(cmds)
            return exe(cmd)
        
        if self.njobs == 0:
            outRet = []
            for c in cmds:
                outRet.append(exe(self.buildCommand(c)))
            return outRet
        
        if id is None:
            id = tempfile.mkstemp(dir=wDir)[1]
        
        outputRet =[]
        for chunk, commands in enumerate( partition(cmds, self.njobs) ):
            outScript = open(os.path.join(wDir, "%s_chunk%d.sh" % (id, chunk)),'w')
            outScript.write("#!/bin/bash\n\n")
            for c in commands:
                outScript.write(c.cmd+"\n")
            outScript.close()
            #Add executeable 
            existing_permissions = stat.S_IMODE(os.stat(outScript.name).st_mode)
            if not os.access(outScript.name, os.X_OK):
                new_permissions = existing_permissions | stat.S_IXUSR
                os.chmod(outScript.name, new_permissions)
                
            submit = Command(outScript.name, \
                            id + "_chunk%d" % chunk, \
                            os.path.join(wDir, id + ("_chunk%d.out" % chunk)), \
                            os.path.join(wDir, id + ("_chunk%d.err" % chunk)))
            cmd = self.buildCommand(submit)
            outputRet.append(exe(cmd))
            
        return outputRet
        
    def checkTemplate(self):
        """
        Checks that my template works okay
        """
        temp.update({"CMD":"test", \
                     "STDOUT":"testo", \
                     "STDERR":"teste", \
                     "JOBNAME":"testn"})
        try:
            w = self.template.substitute(temp)
        except KeyError:
            logging.error("Your submission template is invalid ")
            sys.exit(1)

    def buildCommand(self, cmdSetup):
        """
        substitutes a template with a Command
        """
        return self.template.substitute(cmdSetup.asDict())

def partition(n,m):
    """
    Helper function. splits list n into m partitions
    """
    p = map(lambda x: list(), range(m))
    index = 0
    for item in n:
        p[index].append(item)
        if index < m-1:
            index += 1
        else:
            index = 0
    return filter(lambda x: len(x)>0, p)
