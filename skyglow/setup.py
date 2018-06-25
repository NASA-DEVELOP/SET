from setuptools import setup
setup(
    name = 'skyglow',
    version = '1.0',
    packages = ['skyglow'],
    entry_points = {
        'console_scripts': [
            'skyglow = skyglow.__main__:main'
        ]
    })
