#!/usr/bin/python3
import os,shutil,glob,sys,argparse
import tarfile
from subprocess import Popen, PIPE


def main( a ):
    fo = open('pycrtm.stdo','w')
    fe = open('pycrtm.stde','w')
 
    arch = a.arch
    compilerFlags = selectCompilerFlags(arch)
    # set the required environment variables
    os.environ["FC"] = compilerFlags[arch]['Compiler']  
    os.environ['FCFLAGS']= compilerFlags[arch]['FCFLAGS1']
    installPath = a.install
    crtmRepos = a.rtpath
    scriptDir = os.path.split(os.path.abspath(__file__))[0]
    # if we want to build crtm first.
    if(a.rtinstall):
        # go into crtm git repository directory
        os.chdir( crtmRepos ) 

        print("Configuring/Compiling/Installing CRTM.")
        # configure, comile and install CRTM to the installPath
        configureCompileInstallCrtm( installPath, fo, fe, scriptDir )
        path2CRTM = glob.glob(os.path.join(installPath,'crtm_v*'))

        print("Copying coefficients to {}".format( os.path.join(scriptDir,'crtm_coef_pycrtm') ) )   
        # make the coef directory along with the install location

        # copy coefficients 
        moveCrtmCoefficients( scriptDir  )

        # go back to script directory.
        os.chdir( scriptDir )
    else:
        path2CRTM = installPath
        os.chdir( crtmRepos )
        moveCrtmCoefficients( scriptDir  )
    print("Modifying crtm.cfg")
    modifyOptionsCfg( 'crtm.cfg', scriptDir )
    print("Making python module.")
    # build python module
    # Set compile environment variables.
    os.environ['FCFLAGS'] = os.environ['FCFLAGS'] + compilerFlags[arch]['FCFLAGS2']
    os.environ['FFLAGS'] = os.environ['FCFLAGS']
    os.environ['FC'] = compilerFlags[arch]['Compiler']
    os.environ['ILOC'] = path2CRTM
    os.environ['F2PY_COMPILER'] = compilerFlags[arch]['F2PY_COMPILER']
    os.environ['FORT'] = compilerFlags[arch]['Compiler']
    os.environ['NCINSTALL'] = a.ncpath
    os.environ['H5INSTALL'] = a.h5path
    os.environ['DASHL'] = ' -lcrtm -lgomp -lnetcdf -lnetcdff -lhdf5'
    makeModule(fo, fe, scriptDir)
    os.chdir(scriptDir)
    modifyOptionsCfg( 'crtm.cfg', scriptDir ) 
    fo.close()
    fe.close()
    print("Done!")


def selectCompilerFlags(arch):
    compilerFlags = {}
    if(arch =='gfortran-openmp'):
        compilerFlags['gfortran-openmp']={}
        compilerFlags['gfortran-openmp']['Compiler']='gfortran'
        fullGfortranPath = which('gfortran')
        if(fullGfortranPath ==''): sys.exit("No gfortran found in path.")

        gccBinPath = os.path.split(fullGfortranPath)[0]
        gccPath = os.path.split(gccBinPath)[0]

        # bit to check what gcc version is available, if not > 6. Problem. exit.
        p = Popen(['gfortran','-dumpversion'], stdout = PIPE, stderr = PIPE) 
        p.wait()
        so,se = p.communicate() 
        if ( int(so.decode("utf-8").split('.')[0]) < 6 ):
            sys.exit("F2008 required. gcc >= 6")
        compilerFlags['gfortran-openmp']['FCFLAGS1']="-O3 -fimplicit-none -ffree-form -fno-second-underscore -frecord-marker=4 -funroll-loops -fopenmp -Wall -Wconversion -mieee-fp -fbounds-check -std=f2008"
        compilerFlags['gfortran-openmp']['FCFLAGS2']=""   #mac OS (brew install)
        compilerFlags['gfortran-openmp']['LDFLAGS']="-Wall -g -shared -lnetcdf -lnetcdff -lhdf5"
        compilerFlags['gfortran-openmp']['F2PY_COMPILER']="gnu95"
   
    elif(arch == 'ifort-openmp'):
        compilerFlags['ifort-openmp']={}
        compilerFlags['ifort-openmp']['Compiler']='ifort'
        fullIfortPath = which('ifort')

        if(fullIfortPath == ''): sys.exit("No ifort found.")

        compilerFlags['ifort-openmp']['FCFLAGS1']="-O3 -fp-model source -e08 -free -qopenmp -assume byterecl,realloc_lhs"
        compilerFlags['ifort-openmp']['FCFLAGS2']=" -liomp5 "
        compilerFlags['ifort-openmp']['LDFLAGS']="-Wall -g -shared -liomp5"
        compilerFlags['ifort-openmp']['F2PY_COMPILER']='intelem'
    else:
        sys.exit('Unknown compiler {}.'.format(arch))   
    return compilerFlags
    
def runAndCheckProcess(p, name, fo, fe, scriptDir):
    if(p.returncode>0):
        foname = fo.name
        fename = fe.name

        with open(os.path.join(scriptDir,foname),'r') as foOb:
            for l in foOb.readlines():
                print( l.strip() )

        with open(os.path.join(scriptDir,fename),'r') as feOb:
            for l in feOb.readlines():
                print( l.strip() )

        print("For more information about the install look in {}, and {}".format(fo.name,fe.name) )
        fo.close()
        fe.close()
        sys.exit(name+" failed.")


def configureCompileInstallCrtm( installLocation, fo, fe, scriptDir ):
    # configure as one usually does
    p = Popen(['./configure','--prefix='+installLocation,'--disable-big-endian'],stderr=fe,stdout=fo)
    p.wait()
    runAndCheckProcess(p,"CRTM configure", fo, fe, scriptDir)

    p = Popen(['make','-j'+a.jproc],stderr=fe,stdout=fo,shell=True)
    p.wait()
    runAndCheckProcess(p, "Comipling CRTM", fo, fe, scriptDir)

    # skip make check copying broken in 2.4 alpha release, I think...
    #p = Popen(['make', 'check'],stderr=fe,stdout=fo, shell=True)
    #p.wait()
    #runAndCheckProcess(p,"CRTM check", fo, fe, scriptDir)
    
    p = Popen(['make','install'],stderr=fe,stdout=fo)
    p.wait()        
    runAndCheckProcess(p,"CRTM install", fo, fe, scriptDir)

def moveCrtmCoefficients(installLocation):
    
    if( not os.path.isdir(os.path.join( installLocation, 'crtm_coef_pycrtm' ) ) ):
        os.makedirs( os.path.join( installLocation, 'crtm_coef_pycrtm' ) )
    cwd = os.getcwd()
    p = os.path.join(cwd,'fix','SpcCoeff','Little_Endian')
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','TauCoeff','ODPS','Little_Endian')
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','CloudCoeff','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))
    
    p = os.path.join(cwd,'fix','AerosolCoeff','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','IR_Ice','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','IR_Land','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','IR_Snow','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','IR_Water','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','MW_Water','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','VIS_Ice','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','VIS_Land','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','VIS_Snow','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

    p = os.path.join(cwd,'fix','EmisCoeff','VIS_Water','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), os.path.join(installLocation,'crtm_coef_pycrtm'))

def makeModule(fo, fe, scriptDir):
    # make pycrtm module
    os.chdir(scriptDir)

    p=Popen(['make', 'clean'],stderr=fe,stdout=fo)
    p.wait()
    runAndCheckProcess(p,'pycrtm make clean', fo, fe, scriptDir)
   
    p=Popen(['make'],stderr=fe,stdout=fo)
    p.wait()
    runAndCheckProcess(p,'pycrtm make', fo, fe, scriptDir)

def modifyOptionsCfg( filename, scriptDir ):
    with open( filename+'new' ,'w') as newFile:
        with open(os.path.join(scriptDir, filename), 'r') as oldFile:
            for l in oldFile:
                if('coeffs_dir' in l):
                    newFile.write(l.replace(l,'coeffs_dir = '+os.path.join(os.path.join(scriptDir,'crtm_coef_pycrtm'),'')+os.linesep))
                else:
                    newFile.write(l)
    os.rename(filename+'new', filename)
def which(name):
    found = 0 
    for path in os.getenv("PATH").split(os.path.pathsep):
        full_path = path + os.sep + name
        if os.path.exists(full_path):
            found = 1
            return full_path
    return ''
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = "install pycrtm, optionally crtm if it isn't built already.")
    parser.add_argument('--install',help = 'install path to compiled CRTM install.', required = True, dest='install')
    parser.add_argument('--rtpath',help = 'path to CRTM github clone.', required = True, dest='rtpath')
    parser.add_argument('--ncpath',help = 'path to NETCDF install.', required = True, dest='ncpath')
    parser.add_argument('--h5path',help = 'path to H5 install.', required = True, dest='h5path')
    parser.add_argument('--jproc',help = 'Number of threads to pass to make.', required = True, dest='jproc')
    parser.add_argument('--arch',help = 'compiler/architecture.', required = False, dest='arch', default='gfortran-openmp')
    parser.add_argument('--inplace', help="Switch installer to use rtpath for previously installed crtm.", dest='rtinstall', action='store_false' )
    a = parser.parse_args()
    main(a)
