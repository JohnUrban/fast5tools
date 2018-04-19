import os, tarfile, shutil
from glob import glob
from random import randint
TMP_DIR = ".fast5tools_tmp_dir"
TMP_DIR = TMP_DIR[1:] ## not hdden for now

class FileList(object):
    ## will parse list of information on files, expand FOFNs and Tarballs if necessary, etc
    ## then returns one filename at a time for program to do with as needed.
    def __init__(self, filelist, extension, tar_filenames_only=False, keep_tar_footprint_small=True, downsample=False, random=False, randomseed=False):
        #extension is final 'word' of filename (best including final dot of filename, but not nec) -- e.g. .pdf, .txt, .bam
        #if empty string given, it will then look at all files
        # HOWEVER - with empty string, it will not perceive file.fofn and tarballs as special files to expand
        #ensure type is list
        if isinstance(filelist, list):
                self.filelist = filelist
        elif isinstance(filelist, str):
                self.filelist = [filelist]
        self.extension = tuple(extension) ## can give a single string or list/tuple
        self._tars_detected = False
        self.tar_filenames_only = tar_filenames_only
        self.keep_tar_footprint_small = keep_tar_footprint_small
        self.nfiles = None
##        self.allfiles = None
        self.iterfiles = None
        self.n_tars = 0
        self.TMP_DIR = None
        self.tars = {} ## for small footprint method
        self.files = None
        self.iterfiles = None
        self.most_recent = None
        ## Pre-process
        self.downsample = downsample
        self.random = random
        self.randomseed = randomseed 
        self._extract_files()
        
    def __iter__(self):
        return self

    def next(self):
##        try:
##            return os.path.abspath(self.iterfiles.next())
##        except Exception as e:
##            if self._tars_detected and os.path.exists(TMP_DIR):
##                shutil.rmtree(TMP_DIR)
##            raise StopIteration
        try:
            newfile = self.iterfiles.next()
            if self.keep_tar_footprint_small and newfile.startswith("tar|"):
                self._keep_tar_tmp_dir_small() ## erases last/most_recent file in TMP_DIR
                tar, key, tar_member = newfile.split("|")
                tarkey = "tar|" + key + "|"
                self.tars[tarkey].extract(tar_member, path=self.TMP_DIR)
                newfile = os.path.join(self.TMP_DIR, tar_member)
                newfile = os.path.abspath(newfile)
##                os.remove(newfile)
            else:
                newfile = os.path.abspath(newfile)
            self.most_recent = newfile
            return newfile
        
        except Exception as e:
            if self._tars_detected and os.path.exists(self.TMP_DIR):
                shutil.rmtree(self.TMP_DIR) ## get rid of TMP DIR altogether
            raise StopIteration

    def _initialize_tar_tmp_dir(self):
##        if os.path.isdir(TMP_DIR):
##            shutil.rmtree(TMP_DIR)
##        self._tars_detected = True
##        os.mkdir(TMP_DIR)
        if self.TMP_DIR == None:
            seed = str(randint(10000000000,99999999999))
            self.TMP_DIR = TMP_DIR + "_" + seed
        if os.path.isdir(self.TMP_DIR):
            shutil.rmtree(self.TMP_DIR)
        self._tars_detected = True
        os.mkdir(self.TMP_DIR)

    def _keep_tar_tmp_dir_small(self):
        if self.TMP_DIR == None:
            seed = str(randint(10000000000,99999999999))
            self.TMP_DIR = TMP_DIR + "_" + seed
        if self.most_recent is not None:
            if self.TMP_DIR in self.most_recent:
                os.remove(self.most_recent)

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
        self.n_tars += 1
        if self.tar_filenames_only:
            files = [os.path.join(tarball,fname) for fname in f.getnames() if fname.endswith(self.extension)]
        else: ## will be using files in tarball
            if self.keep_tar_footprint_small:
                tarkey = "tar|"+str(self.n_tars)+"|"
                self.tars[tarkey] = f
                files = [tarkey+fname for fname in f.getnames() if fname.endswith(self.extension)]
                ## purposely do not close tarfile
            else:
                f.extractall(path=TMP_DIR) ## in situations where tarball includes many big non-file files, this may not be best way.
                files = [os.path.join(TMP_DIR, fname) for fname in f.getnames() if fname.endswith(self.extension)]
        if not self.keep_tar_footprint_small: ## need to keep tarchive open otherwise
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
        if self.downsample:
            if self.random:
                self.randomseed = self.randomseed if self.randomseed else randint(0,1000000)
            self.down_sample_iter_files(n=self.downsample, random=self.random, randomseed=self.randomseed, sort=True)

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


    def down_sample_iter_files(self, n=1, random=False, randomseed=False, sort=True):
        if n >= self.nfiles or n <= 0 or n is None or n is False: ## no downsampling possible, return all
            self.iterfiles =  iter(self.files)
        else:
            self.downsampled_files = self.files[:]
            if random:
                if randomseed:
                    seed(randomseed)
                shuffle(self.downsampled_files)
            self.downsampled_files = self.downsampled_files[:n]
            if sort:
                ## Just trying to return the sampled, potentially random list as sorted
                self.downsampled_files = sorted(self.downsampled_files)
            #print self.downsampled_files
            self.iterfiles =  iter(self.downsampled_files)

    def reset_iter_files(self):
        self.iterfiles =  iter(self.files)
