import os, tarfile
from glob import glob
TMP_DIR = ".fast5tools_tmp_dir"
TMP_DIR = TMP_DIR[1:]
class FileList(object):
    ## will parse list of information on files, expand FOFNs and Tarballs if necessary, etc
    ## then returns one filename at a time for program to do with as needed.
    def __init__(self, filelist, extension, tar_filenames_only=False):
        #extension is final 'word' of filename (best including final dot of filename, but not nec) -- e.g. .pdf, .txt, .bam
        #if empty string given, it will then look at all files
        # HOWEVER - with empty string, it will not perceive file.fofn and tarballs as special files to expand
        #ensure type is list
        if isinstance(filelist, list):
                self.filelist = filelist
        elif isinstance(filelist, str):
                self.filelist = [filelist]
        self.extension = extension
        self._tars_detected = False
        self.tar_filenames_only = tar_filenames_only
        self.nfiles = None
##        self.allfiles = None
        self.files = None
        self.iterfiles = None
        self._extract_files()
        
    def __iter__(self):
        return self

    def next(self):
        try:
            return os.path.abspath(self.iterfiles.next())
        except Exception as e:
            if self._tars_detected and os.path.exists(TMP_DIR):
                shutil.rmtree(TMP_DIR)
            raise StopIteration

    def _initialize_tar_tmp_dir(self):
        if os.path.isdir(TMP_DIR):
            shutil.rmtree(TMP_DIR)
        self._tars_detected = True
        os.mkdir(TMP_DIR)

    def _expand_fofn(self, fofn):
        f = open(fofn,'r')
        lines = [line.strip() for line in f.readlines()]
        f.close()
        files = []
        for fname in lines:
            if fname.endswith(self.extension):
                files.append(fname)
            elif os.path.isdir(fname):
                files += self._expand_dir(fname)
            elif tarfile.is_tarfile(fname):
                files += self._expand_tar(fname)
            else: ## line in FOFN is ignored
                pass
        return files

    def _expand_dir(self, d):
        pattern = d + '/' + '*'+self.extension
        files = glob(pattern)
        return files

    def _expand_tar(self, tarball):
        if not self._tars_detected and not self.tar_filenames_only:
            self._initialize_tar_tmp_dir()
        f = tarfile.open(tarball)
        if self.tar_filenames_only:
            files = [os.path.join(tarball,fname) for fname in f.getnames() if fname.endswith(self.extension)]
        else: ## will be using files in tarball
            f.extractall(path=TMP_DIR) ## in situations where tarball includes many big non-file files, this may not be best way.
            files = [os.path.join(TMP_DIR, fname) for fname in f.getnames() if fname.endswith(self.extension)]
        f.close()
        return files

    def _extract_files(self):
        self.files = []
        for e in self.filelist:
            if e.endswith(self.extension):
                self.files.append(e)
            elif e.endswith(".fofn"):
                self.files += self._expand_fofn(e)
            elif os.path.isdir(e):
                self.files += self._expand_dir(e)
            elif tarfile.is_tarfile(e):
                self.files += self._expand_tar(e)
            else: ## all else currently ignored 
                pass
        self.nfiles = len(self.files)
        self.iterfiles = iter(self.files)

    def get_file_list(self):
        return self.files
    
    def __len__(self):
        return self.nfiles

    def __eq__(self, other):
        return set(self.files) == set(other.files)

    def get_filenames(self):
        return self.files

    def get_basenames(self):
        return [os.path.basename(e) for e in self.files]

    def get_dirnames(self):
        return [os.path.dirname(e) for e in self.files]
