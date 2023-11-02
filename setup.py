import os, sys, tarfile, configparser, ssl, shutil, socket
from contextlib import closing
import urllib.request
from skbuild import setup
def main():
    #Completely remove previous _skbuild, because cache will remember previous interation and ignore you if you change something.
    try:shutil.rmtree('_skbuild')
    except:pass 
    
    #path of this file.
    scriptDir = os.path.split(os.path.abspath(__file__))[0]
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
        packages=['pycrtm'],
        py_modules=['crtm_io', 'pyCRTM'],
        package_data={'pyCRTM':['pyCRTM/setup.txt']})
    
    os.remove('MANIFEST.in')
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
        os.makedirs( coefDest )
    else:
        shutil.rmtree( coefDest )
        os.makedirs( coefDest )

    topdir = os.listdir(coefDir)
    # if its more than 1 item, this means we've hit a "fix" directory
    # otherwise treat it like an extracted coef tarball
    if(len(topdir)>1):
        for t in topdir:
            if( os.path.isfile( os.path.join(coefDir,t) ) ):
                os.symlink( os.path.join(coefDir,t), os.path.join(coefDest,t))
    else:
        td = topdir[0]
        searchPath = os.path.join(coefDir,td)
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

        for l in toLink:
            os.symlink(l, os.path.join(coefDest,os.path.split(l)[1]))
if __name__ == "__main__":
    main()

