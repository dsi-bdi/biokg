from zipfile import ZipFile, ZIP_DEFLATED
from os import walk
from os.path import join, basename

data_root = 'data/'


def package_biokg():
    biokg = ZipFile(join(data_root, 'biokg.zip'), 'w', compression=ZIP_DEFLATED, compresslevel=6)
    biokg.write('LICENSE')
    biokg.write('README.md')

    for folderName, subfolders, filenames in walk(join(data_root, 'biokg')):
        for filename in filenames:
            #create complete filepath of file in directory
            filePath = join(folderName, filename)
            # Add file to zip
            biokg.write(filePath, basename(filePath))
    biokg.close()


def package_benchmarks():
    bmark = ZipFile(join(data_root, 'benchmarks.zip'), 'w', compression=ZIP_DEFLATED, compresslevel=6)
    bmark.write('LICENSE')

    for folderName, subfolders, filenames in walk(join(data_root, 'output', 'benchmarks')):
        for filename in filenames:
            #create complete filepath of file in directory
            filePath = join(folderName, filename)
            # Add file to zip
            bmark.write(filePath, basename(filePath))
    bmark.close()


def package_all():
    package_biokg()
    package_benchmarks()


if __name__ == '__main__':
    package_all()