import os, sys, configparser, shutil, time
from skbuild import setup
def main():
    #Completely remove previous _skbuild, because cache will remember previous interation and ignore you if you change something.
    try:shutil.rmtree('_skbuild')
    except:pass 
    
    #path of this file.
    scriptDir = os.path.split(os.path.abspath(__file__))[0]
    #add .F90 symlink to keep history of pycrtm.f90 
    if( not os.path.exists( os.path.join(scriptDir,'pycrtm.F90'))):
        os.symlink(os.path.join(scriptDir,'pycrtm.f90'),os.path.join(scriptDir,'pycrtm.F90'))
    #read configuration
    coef_path, coef_dest, crtm_install, link_coef = readSetup('setup.cfg',scriptDir)
    if(link_coef):
        linkCoef(coef_path, coef_dest)
    os.environ['CRTM_INSTALL'] = crtm_install
    shutil.copy(os.path.join(scriptDir,'setup.cfg'),os.path.join(scriptDir,'pyCRTM','pycrtm_setup.txt'))

    f = open(os.path.join(scriptDir,'MANIFEST.in'),'w')
    f.write('include pyCRTM/pycrtm_setup.txt crtm_io.py')
    f.close()

    requires=['numpy']
    setup(
        name="pyCRTM_JCSDA",
        version='2.0.1',
        description='Python wrapper for the CRTM.',
        author='Bryan Karpowicz',
        requires=requires,
        include_package_data=True,
        packages=['pycrtm_'],
        py_modules=['crtm_io', 'pyCRTM'],
        package_data={'pyCRTM':['pyCRTM/setup.txt']})
    
    os.remove('MANIFEST.in')
    #os.remove('pycrtm.F90')
def readSetup(setup_file, scriptDir):
    cfg = configparser.ConfigParser()
    cfg.read( os.path.join(scriptDir,'setup.cfg') )
    crtm_install = cfg['Setup']['crtm_install']
    coef_path = cfg['Coefficients']['source_path']
    coef_dest = cfg['Coefficients']['path_used']
    link_coef = cfg['Setup']['link_from_source_to_path_used']
    return coef_path, coef_dest, crtm_install, link_coef

def linkCoef(coefDir,coefDest):
    print("Linking Coefficients.") 
    cwd = os.getcwd()

    if( not os.path.isdir( coefDest )  ):
        print("Creating directory:{}".format(coefDest))
        os.makedirs( coefDest )
    else:
        print("Cleaning directory:{}".format(coefDest))
        shutil.rmtree( coefDest )
        print("Creating directory:{}".format(coefDest))
        os.makedirs( coefDest )
    print("Linking Coefficients in {} to {}".format(coefDir,coefDest))
    topdir = os.listdir(coefDir)
    # use os.walk to get files regardless if they're in an emc "fix" style directory, or crtm coeff tarball.
    searchPath = coefDir 
    toLink = []
    filesPresent = []
    for root,dirs,filez in os.walk(searchPath):
        for name in filez:
            srcPath = os.path.join(root,name)
            if(not ('ODAS' in srcPath or 'Big_Endian' in srcPath)):
                curF = os.path.split(srcPath)[1]
                if (curF not in filesPresent):
                    toLink.append(srcPath)
                    filesPresent.append(curF)
    coefCnt = 0
    for l in toLink:
        os.symlink(l, os.path.join(coefDest,os.path.split(l)[1]))
        coefCnt+=1
    if coefCnt == 0:
        print("Warning! Linked zero coefficients!")
        time.sleep(30)
    else:
        print("Linked {:d} Coefficients.".format(coefCnt)) 
    time.sleep(10)
if __name__ == "__main__":
    main()

