import os, sys, configparser, shutil, time, glob
import setuptools
import platform
from skbuild import setup

def is_windows():
    return platform.system() == "Windows"

def copy_file_windows(src, dst):
    """Windows下复制文件而不是创建符号链接"""
    try:
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        print(f"复制文件失败 {src} -> {dst}: {e}")
        return False

def main():
    # 完全移除之前的构建目录
    try:
        shutil.rmtree('_skbuild')
    except:
        pass 
    
    # 当前文件路径
    scriptDir = os.path.split(os.path.abspath(__file__))[0]

    # 读取配置
    coef_path, coef_dest, crtm_install, link_coef = readSetup('setup.cfg', scriptDir)
    
    # Windows特殊处理
    if is_windows():
        print("检测到Windows系统，应用特殊处理...")
        # 使用用户目录避免权限问题
        user_dir = os.path.expanduser("~")
        coef_dest = os.path.join(user_dir, "pycrtm_coefficients")
        crtm_install = coef_dest
        
    if link_coef:
        linkCoef(coef_path, coef_dest)
        
    os.environ['CRTM_INSTALL'] = crtm_install
    
    shutil.copy(os.path.join(scriptDir,'setup.cfg'), os.path.join(scriptDir,'pycrtm_setup.txt'))

    requires=['numpy']
    setup(
        name="pycrtm_jcsda",
        version='2.0.1',
        description='Python wrapper for the CRTM.',
        author='Bryan Karpowicz',
        requires=requires,
        include_package_data=True,
        packages=['pycrtm_'],
        py_modules=['crtm_io', 'pyCRTM'],
        zip_safe=False,
    )

def readSetup(setup_file, scriptDir):
    cfg = configparser.ConfigParser()
    cfg.read(os.path.join(scriptDir, 'setup.cfg'))
    crtm_install = cfg['Setup']['crtm_install']
    coef_path = cfg['Coefficients']['source_path']
    coef_dest = cfg['Coefficients']['path_used']
    link_coef = cfg['Setup']['link_from_source_to_path_used'] == 'True'
    return coef_path, coef_dest, crtm_install, link_coef

def linkCoef(coefDir, coefDest):
    print("处理系数文件...") 
    
    # Windows特殊处理
    if is_windows():
        return linkCoef_windows(coefDir, coefDest)
    else:
        return linkCoef_unix(coefDir, coefDest)

def linkCoef_windows(coefDir, coefDest):
    """Windows版本的系数文件处理"""
    print("Windows系统: 使用文件复制代替符号链接")
    
    if not os.path.isdir(coefDest):
        print(f"创建目录: {coefDest}")
        os.makedirs(coefDest)
    else:
        print(f"清理目录: {coefDest}")
        # 只删除文件，不删除目录
        for filename in os.listdir(coefDest):
            file_path = os.path.join(coefDest, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(f"删除文件失败 {file_path}: {e}")

    print(f"复制系数文件从 {coefDir} 到 {coefDest}")
    
    # 搜索所有系数文件
    searchPath = coefDir 
    toCopy = []
    filesPresent = []
    
    for root, dirs, filez in os.walk(searchPath):
        for name in filez:
            srcPath = os.path.join(root, name)
            # 跳过不需要的目录
            if not ('ODAS' in srcPath or 'Big_Endian' in srcPath):
                curF = os.path.split(srcPath)[1]
                if curF not in filesPresent:
                    toCopy.append(srcPath)
                    filesPresent.append(curF)
    
    # 复制文件
    copyCnt = 0
    for src in toCopy:
        dst = os.path.join(coefDest, os.path.split(src)[1])
        if copy_file_windows(src, dst):
            copyCnt += 1
    
    if copyCnt == 0:
        print("警告: 没有复制任何系数文件!")
        # 在Windows上不暂停，避免阻塞
    else:
        print(f"复制了 {copyCnt} 个系数文件")
    
    # 检查缺失的ODAS系数
    check_missing_coefficients(coefDir, coefDest)
    
    return copyCnt

def linkCoef_unix(coefDir, coefDest):
    """Unix版本的系数文件处理（原逻辑）"""
    print("Linking Coefficients.") 
    
    if not os.path.isdir(coefDest):
        print(f"Creating directory: {coefDest}")
        os.makedirs(coefDest)
    else:
        print(f"Cleaning directory: {coefDest}")
        shutil.rmtree(coefDest)
        print(f"Creating directory: {coefDest}")
        os.makedirs(coefDest)
        
    print(f"Linking Coefficients in {coefDir} to {coefDest}")
    
    # 原版的符号链接逻辑
    toLink = []
    filesPresent = []
    for root, dirs, filez in os.walk(coefDir):
        for name in filez:
            srcPath = os.path.join(root, name)
            if not ('ODAS' in srcPath or 'Big_Endian' in srcPath):
                curF = os.path.split(srcPath)[1]
                if curF not in filesPresent:
                    toLink.append(srcPath)
                    filesPresent.append(curF)
    
    linkCnt = 0
    for src in toLink:
        dst = os.path.join(coefDest, os.path.split(src)[1])
        os.symlink(src, dst)
        linkCnt += 1
        
    if linkCnt == 0:
        print("Warning! Linked zero coefficients!")
        time.sleep(30)
    else:
        print(f"Linked {linkCnt} Coefficients.")
    
    check_missing_coefficients(coefDir, coefDest)
    
    return linkCnt

def check_missing_coefficients(coefDir, coefDest):
    """检查缺失的系数文件"""
    print("Checking for missing ODAS only coefficients.") 
    
    SpcCoeffNc = glob.glob(os.path.join(coefDest, '*SpcCoeff*.nc'))
    SpcCoeffBin = glob.glob(os.path.join(coefDest, '*SpcCoeff*.bin'))
    missingTau = []
    
    for c in SpcCoeffNc:
        if not os.path.exists(c.replace('SpcCoeff', 'TauCoeff')):
            cc = os.path.split(c.replace('SpcCoeff', 'TauCoeff'))[-1]
            missingTau.append(cc)
            
    for c in SpcCoeffBin:
        if not os.path.exists(c.replace('SpcCoeff', 'TauCoeff')):
            cc = os.path.split(c.replace('SpcCoeff', 'TauCoeff'))[-1]
            missingTau.append(cc)
    
    # 处理缺失的文件
    toProcess = []
    filesPresent = []
    
    for root, dirs, filez in os.walk(coefDir):
        for name in filez:
            srcPath = os.path.join(root, name)
            if not ('ODPS' in srcPath or 'Big_Endian' in srcPath):
                curF = os.path.split(srcPath)[1]
                if curF in missingTau and curF not in filesPresent:
                    toProcess.append(srcPath)
                    filesPresent.append(curF)
    
    processedCnt = 0
    for src in toProcess:
        dst = os.path.join(coefDest, os.path.split(src)[1])
        if is_windows():
            if copy_file_windows(src, dst):
                processedCnt += 1
        else:
            os.symlink(src, dst)
            processedCnt += 1
    
    print(f"Processed {processedCnt} additional coefficients.") 
    
    if not is_windows():
        time.sleep(10)

if __name__ == "__main__":
    main()
