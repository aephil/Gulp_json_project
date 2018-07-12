from setuptools import setup

scripts_phononwebsite = ['scripts/read_gulp.py',]
packages_phononwebsite = ['Gulp_json_project']

if __name__ == '__main__':
    setup(name='Gulp_json_project',
          version='0.1',
          description='Read Gulp simulator phonon dispersions to visualize on the phonon website.',
          author='Tushar Bairagi',
          author_email='tusharsanjaybairagi@gmail.com',
          scripts=scripts_phononwebsite,
          packages=packages_phononwebsite
          )


