# -*- coding: utf-8 -*-

import ssl
import hashlib
import os
import tarfile
import urllib.request
import shutil
from tempfile import gettempdir
from timeit import default_timer as timer
from os.path import join, basename, isfile, isdir, dirname
from os import mkdir
import requests
from .extras import *


def download_file_with_cert(url, local_path, checksum=None, cert=None):
    """
    download a file with certificate to disk with ability to validate file checksum

    :param url: string represents file full web url
    :param local_path: string represents full local path
    :param checksum: string represents the checksum of the file
    :param cert: ssl.context context certificate
    :return: string represents local file full path

    This function will download a file to disk from the specified url with checking certificate
    to the specified local file path.
    """

    if cert is None:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    else:
        ctx = cert

    tmp_dir = gettempdir()
    file_name = url.split('/')[-1]
    tmp_file = os.path.join(tmp_dir, file_name)

    with urllib.request.urlopen(url, context=ctx) as u, open(tmp_file, mode='wb') as f:
        f.write(u.read())

    # check the checksum if provided
    if not (checksum is None):
        # compute the checksum of the file
        downloaded_file_checksum = get_file_md5(tmp_file)
        if downloaded_file_checksum != checksum:
            raise ValueError("invalid file checksum [%s] for file: %s" % (downloaded_file_checksum, url))

    # move tmp file to desired file local path
    shutil.move(tmp_file, local_path)
    return local_path


def download_file(url, local_path, checksum=None):
    """ download a file to disk with ability to validate file checksum

    Parameters
    ----------
    url : str
        represents file full web url
    local_path : str
        represents full local path
    checksum : str
        represents the checksum of the file
    Returns
    -------
    str
        local path to downloaded file
    """

    # create a temporary file to download to
    tmp_dir = gettempdir()
    file_name = url.split('/')[-1]
    tmp_file = os.path.join(tmp_dir, file_name)

    header = [('User-agent', 'Mozilla/5.0 (X11; Fedora; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/75.0.3770.100 Safari/537.36')]
    opener = urllib.request.build_opener()
    opener.addheaders = header
    urllib.request.install_opener(opener)
    urllib.request.urlretrieve(url, tmp_file)

    # check the checksum if provided
    if not (checksum is None):
        # compute the checksum of the file
        downloaded_file_checksum = get_file_md5(tmp_file)
        if downloaded_file_checksum != checksum:
            raise ValueError("invalid file checksum [%s] for file: %s" % (downloaded_file_checksum, url))

    # move tmp file to desired file local path
    shutil.move(tmp_file, local_path)
    return local_path


def download_file_with_auth(url, local_path, username, password, checksum=None):
    """ download a file to disk using the given uername and password 
    with ability to validate file checksum

    Parameters
    ----------
    url : str
        represents file full web url
    local_path : str
        represents full local path
    username : str
        the username to use for authentication
    password : str
        the password to use for authentication
    checksum : str
        represents the checksum of the file
    Returns
    -------
    str
        local path to downloaded file
    """

    # create a temporary file to download to
    tmp_dir = gettempdir()
    file_name = url.split('/')[-1]
    tmp_file = os.path.join(tmp_dir, file_name)
    
    with requests.get(url, auth=(username, password)) as r:
        if r.status_code == 200:
            with open(tmp_file, 'wb') as out:
                for bits in r.iter_content():
                    out.write(bits)
        else:
            raise ValueError(f'Error downloading {url} status: {r.status_code}')
    # check the checksum if provided
    if not (checksum is None):
        # compute the checksum of the file
        downloaded_file_checksum = get_file_md5(tmp_file)
        if downloaded_file_checksum != checksum:
            raise ValueError("invalid file checksum [%s] for file: %s" % (downloaded_file_checksum, url))

    # move tmp file to desired file local path
    shutil.move(tmp_file, local_path)
    return local_path


def get_file_md5(file_path):
    """ Get the md5 hash of a file

    Parameters
    ----------
    file_path : str
        represents file path

    Returns
    -------
    str
        The file MD5 hash
    """

    hash_md5 = hashlib.md5()
    with open(file_path, mode='rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def verify_md5(folder, dataset_file):
    """ Verify if the MD5 of a file match the MD5 included in a checksum file.
    The function checks if the computed MD5 of :code:`dataset_file` included in :code:`folder`
    matches the MD5 included in the first row of the file :code:`CHECKSUMS` stored in the same location.

    Parameters
    ----------
    folder : str
        the folder of :code:`dataset_file` and its checksum file.
    dataset_file : str
        the name of the file to check.
    Returns
    -------
    bool
        True if MD5s match, False otherwise.
    """
    actual_md5 = get_file_md5(os.path.join(folder, dataset_file))
    with open(os.path.join(folder, 'CHECKSUMS')) as f_check:
        expected_md5 = f_check.readline().split()[0]
    return True if actual_md5 == expected_md5 else False


def extract_tgz_file(file_path, extraction_dir):
    """ Extract compressed file into a directory

    Parameters
    ----------
    file_path : str
        represents the compressed file path
    extraction_dir : str
        represents the extraction directory
    """

    fd = tarfile.open(file_path)
    fd.extractall(extraction_dir)
    fd.close()


def download_file_md5_check(download_url, filepath, username=None, password=None, bypass_md5=False):
    """ Download file if not existing and have valid md5

    Parameters
    ----------
    download_url : str
        file download url
    filepath : str
        file full local path
    """
    filename = basename(filepath)
    md5_dp = join(dirname(filepath), "checksum")
    md5_filepath = join(md5_dp, filename + ".md5")

    require_download = False
    #File exists with no md5 and bypass_md5 is true
    #Assumes file is valid and creates md5 for file
    if bypass_md5 and isfile(filepath) and not isfile(md5_filepath):
        file_computed_md5 = get_file_md5(filepath)
        mkdir(md5_dp) if not isdir(md5_dp) else None
        md5_fd = open(md5_filepath, "w")
        md5_fd.write(file_computed_md5)
        md5_fd.close()
        print(hsh_sym + "md5 hash saved (%s)." % file_computed_md5, flush=True)
    elif isfile(filepath) and isfile(md5_filepath):
        file_computed_md5 = get_file_md5(filepath)
        file_saved_md5_fd = open(md5_filepath, "r")
        file_saved_md5 = file_saved_md5_fd.read().strip()
        file_saved_md5_fd.close()
        print(inf_sym + "file (%-40s) exists." % filename, end="", flush=True)
        if file_computed_md5 == file_saved_md5:
            print(hsh_sym + " MD5 hash (%s) is correct%s. No download is required." % (file_saved_md5, done_sym))
        else:
            require_download = True
            print(hsh_sym + " MD5 hash (%s) is inconsistent%s. Re-downloading ..." % (file_computed_md5, fail_sym))
    else:
        require_download = True

    if require_download:
        print(dwn_sym + "downloading file (%-40s) ..." % filename, end="", flush=True)
        start = timer()
        if username is None or password is None:
            download_file(download_url, filepath)
        else:
            download_file_with_auth(download_url, filepath, username, password)
        download_time = timer() - start
        print(done_sym + " %1.2f Seconds." % download_time, end="", flush=True)
        file_computed_md5 = get_file_md5(filepath)
        mkdir(md5_dp) if not isdir(md5_dp) else None
        md5_fd = open(md5_filepath, "w")
        md5_fd.write(file_computed_md5)
        md5_fd.close()
        print(hsh_sym + "md5 hash saved (%s)." % file_computed_md5, flush=True)


def export_file_md5(filepath):
    """ Save the file's md5 hash to file

    Parameters
    ----------
    filepath : str
        file full path

    Returns
    -------
    bool
        True if the file has valid md5 and false if not
    """
    computed_md5 = get_file_md5(filepath)
    checksum_dp = join(dirname(filepath), "checksum")
    md5_fp = join(checksum_dp, basename(filepath) + ".md5")
    os.makedirs(checksum_dp) if not os.path.isdir(checksum_dp) else None
    md5_fd = open(md5_fp, "w")
    md5_fd.write(computed_md5)
    md5_fd.close()


def file_has_valid_md5(filepath):
    """ Check that the file has valid md5

    Parameters
    ----------
    filepath : str
        file full path

    Returns
    -------
    bool
        True if the file has valid md5 and false if not
    """
    md5_fp = join(dirname(filepath), "checksum", basename(filepath) + ".md5")
    if isfile(filepath) and isfile(md5_fp):
        saved_md5 = open(md5_fp).read().strip()
        computed_md5 = get_file_md5(filepath)
        if saved_md5 == computed_md5:
            return True
        else:
            return False
    else:
        return False
