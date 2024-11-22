#!/usr/bin/env python

##
#  Retrieves the SCEC project location from CARC. They are
#  very big
#

import getopt
import sys
import subprocess
import os

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen

model = "CS248"

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url,fname):
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True

def main():

    # Set our variable defaults.
    fname = ""
    path = ""
    mdir = ""

    # Get the dataset.
    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)

    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()
        if (variable == 'model_data_path') :
            path = val + '/' + model
            continue
        if (variable == 'model_dir') :
            mdir = "./"+val
            continue

        continue
    if path == "" :
        print("ERROR: failed to find variables from config file")
        sys.exit(1)

    fp.close()

    print("\nDownloading model dataset\n")


#check if cvm-large-dataset/cs248 exists or not
#yes,  link it over
    volume_top_dir=os.environ.get('CVM_VOLUME_TOP_DIR')
    if volume_top_dir != None :
        dpath=volume_top_dir+"/"+mdir
        if os.path.isdir(dpath) :
            subprocess.check_call(["ln", "-s", dpath, "."])
            print("\nLinked!")
            return;

#no, then download only if needs too.
    if not os.path.isdir(mdir) :
        subprocess.check_call(["mkdir", "-p", mdir])

    fname=mdir+"/density.dat"
    tarfile=mdir+"/density.tar.gz"
    if not os.path.isfile(fname) :
      print("download ", tarfile)
      url = path + "/" + tarfile 
      download_urlfile(url,tarfile)
      subprocess.check_call(["tar", "-zxvf", tarfile])

    fname=mdir+"/vp.dat"
    tarfile=mdir+"/vp.tar.gz"
    if not os.path.isfile(fname) :
      print("download ", tarfile)
      url = path + "/" + tarfile
      download_urlfile(url,tarfile)
      subprocess.check_call(["tar", "-zxvf", tarfile])

    fname=mdir+"/vs.dat"
    tarfile=mdir+"/vs.tar.gz"
    if not os.path.isfile(fname) :
      print("download ", tarfile)
      url = path + "/" + tarfile
      download_urlfile(url,tarfile)
      subprocess.check_call(["tar", "-zxvf", tarfile])

    print("\nDone!")

if __name__ == "__main__":
    main()
