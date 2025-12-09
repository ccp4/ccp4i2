import sys
import zipfile

def extract_results_zip(fn):
    if not zipfile.is_zipfile(fn):
        return
    print("OK, a zip file")
    with zipfile.ZipFile(fn) as myzip:
       infolist = myzip.infolist()
       for finfo in infolist:
           if finfo.filename.endswith("_besttls.mtz"):
               print(finfo.filename)
               myzip.extract(finfo)
           if finfo.filename.endswith("_final.mtz"):
               print(finfo.filename)
               myzip.extract(finfo)
           if finfo.filename.endswith("_besttls.pdb"):
               print(finfo.filename)
               myzip.extract(finfo)
           if finfo.filename.endswith("_final.pdb"):
               print(finfo.filename)
               myzip.extract(finfo)

if __name__ == "__main__":
    fn = sys.argv[1]
    extract_results_zip(fn)
