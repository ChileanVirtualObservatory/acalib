import sys
import time
import urllib
import os

from asydopy import *

default_url = "http://www.csrg.cl/~maray/splatalogue.csv"
default_csv_name = "lines2.csv"
default_db_name = "ASYDO"
csv = False
URI = ""
log = sys.stdout

database = db.lineDB(default_db_name)

def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    log.write("\r\t...%d%%, %d MB, %d KB/s, %d seconds passed" % (percent, progress_size / (1024 * 1024), speed, duration))
    log.flush()

def helper(log):
    log.write("\n")
    log.write("usage:   dbCreator.py | [options] [path] | [options] [path] [secondary] [range]| \n")
    log.write("Options:\n")
    log.write("\t -C\t: Use .CSV as input file, [path] is path to the desired file \n")
    log.write("\t -T\t: Use SLAP Service as input, [path] is URL of service \n")
    log.write("Secundary:\n")
    log.write("\t -R\t: Only used wit -T, allows for custom wavelengh range query. [range] is in the form <minimum:maximum>\n")
    log.write("\n")

complete = False

if (len(sys.argv)==1):

    log.write("No source was defined for the Line Database, so the Default Line Database will be imported: \n\tDownloading Default Line CSV \n")
    urllib.urlretrieve(default_url,default_csv_name , reporthook)
    URI = default_csv_name
    log.write("\n\tDefault Line CSV was Downloaded\n")
    csv = True

elif (len(sys.argv)>2):
    URI = sys.argv[2]

    if (sys.argv[1]=="-C"):
        log.write("Using %s as source for Line Database\n" % os.path.basename(URI))
        csv = True

    elif sys.argv[1]=="-T":
        w_range = [88,720]
        message = "No range specified, using Default range (88 Ghz to 720 Ghz)\n"
        if len(sys.argv) == 5 and sys.argv[3] == "-R":
            w_range = sys.argv[4].split(":")
            message = "The specified range is ("+w_range[0]+" Mhz to " + w_range[1]+" Mhz)\n"

        log.write(message)
        database.deleteDB()
        database.VOGetLines(log,URI,w_range)
        database.loadVoTable("./votables/customVOTable.xml",{3:"FREQ",4:"SPECIES",5:"CHEM_NAME",7:"INTENSITY",11:"EL"})
        complete = True


    else:
        helper(log)
else:
    helper(log)


if csv:
    log.write("Importing CSV (%s) to SQL Database\n" % os.path.basename(URI))
    database.createDBFromCSV(URI,log)
    complete = True

if complete:
    log.write("Database creation is now complete, you can now use ASYDO. Have a nice day.\n")
