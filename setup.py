from setuptools import setup, find_packages


setup(
    name='PacBE',
    package_dir={'': 'src'},      
    packages=find_packages(where='src'),
    include_package_data=True,
	py_modules=['PacBe_demultiplex', 
		'PacBe_organize_data'],
    entry_points={
        'console_scripts': [
            'PacBe_demultiplex = PacBe_demultiplex:main',
            'PacBe_organize_data = PacBe_organize_data:main'
        ],
    },
    
)
