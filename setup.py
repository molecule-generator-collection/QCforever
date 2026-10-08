"""QCforever is a wrapper of Gaussian.    
    
QCforever (https://doi.org/10.1021/acs.jcim.2c00812) is a wrapper of Gaussian (https://gaussian.com).    
To compute obsevable properties of a molecule through quantum chemical computation (QC), Multi step computation is demanded.
QCforever automates this process and calculates multiple physical properties of molecules simultaneously.                       
""" 
from setuptools import setup
                                
    
DOCLINES = (__doc__ or '').split('\n')
INSTALL_REQUIRES = [                      
    'numpy>=1.22.0',    
    'rdkit>=2023.3.3',
    'basis-set-exchange>=0.12',
    'bayesian-optimization==1.4.3',
    'psutil']
INSTALL_REQUIRES.append('PyYAML>=6')
PACKAGES = [                    
    'qcforever',
    'qcforever.gaussian_run',
    'qcforever.gamess_run',
    'qcforever.util',
    'qcforever.laqa_fafoom',
    'qcforever.conformer_search',
    'qcforever.conformer_search.model_workers']
PACKAGE_DATA = {
    'qcforever.conformer_search.model_workers': ['requirements/*.txt', 'model_sources.json'],
    'qcforever.conformer_search': ['defaults.yaml'],
    'qcforever': ['gaussian_run/*.json', 'gaussian_run/*.csv'],
}
CLASSIFIERS = [           
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.11',    
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
    version="2.5.0",                                                     
    download_url="https://github.com/molecule-generator-collection/QCforever",
    python_requires=">=3.11",                                                      
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,                    
    package_data=PACKAGE_DATA,
    entry_points={'console_scripts': [
        'qcforever-model-worker=qcforever.conformer_search.model_workers.run_model:main',
        'install-conformer-models=qcforever.conformer_search.model_workers.install_models:main',
    ]},
    classifiers=CLASSIFIERS
)  
