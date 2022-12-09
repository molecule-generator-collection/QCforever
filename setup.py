"""QCforever is a wrapper of Gaussian.    
    
QCforever (https://doi.org/10.1021/acs.jcim.2c00812) is a wrapper of Gaussian (https://gaussian.com).    
To compute obsevable properties of a molecule through quantum chemical computation (QC), Multi step computation is demanded.
QCforever automates this process and calculates multiple physical properties of molecules simultaneously.                       
""" 
from setuptools import setup
                                
    
DOCLINES = (__doc__ or '').split('\n')
INSTALL_REQUIRES = [                      
    'numpy~=1.19.2',    
    'rdkit-pypi==2021.03.5']
PACKAGES = [                    
    'qcforever',
    'qcforever.gaussian_run']
PACKAGE_DATA = {
    'qcforever': ['gaussian_run/*.json'],
}
CLASSIFIERS = [           
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',    
    "License :: OSI Approved :: MIT License",   
    "Operating System :: OS Independent"]        
                                             
setup(
    name="qcforever",
    author="Masato Sumita",
    author_email="masato.sumita@riken.jp",
    maintainer="Masato Sumita",               
    maintainer_email="masato.sumita@riken.jp",
    description=DOCLINES[0],                      
    long_description='\n'.join(DOCLINES[2:]),
    long_description_content_type="text/markdown",
    license="MIT LIcense",                            
    url="https://github.com/molecule-generator-collection/QCforever",
    version="1.0.0",                                                     
    download_url="https://github.com/molecule-generator-collection/QCforever",
    python_requires=">=3.7",                                                      
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,                    
    package_data=PACKAGE_DATA,
    classifiers=CLASSIFIERS
)  