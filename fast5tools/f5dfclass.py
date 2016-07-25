import numpy as np
import pandas

class Fast5DataFrame(object):
    def __init__(self, df):
        '''df is a pandas dataframe object'''
        self.df = df
        self.has2d = self.df['has2d'] == 1
        self.hascomp = self.df['hascomp'] == 1
        self.no2d = self.df['has2d'] == 0
        self.temponly = self.df['hascomp'] == 0
        self.n_molecules = None
        self.n_temp = None
        self.n_temp_only = None
        self.n_comp = None
        self.n_no_2d = None
        self.n_comp_has_2d = None
        self.n_comp_no_2d = None
        self.n_2d = None
        self.sum_molecule_lengths = None
        self.nx_mol = None
        self.lx_mol = None
        self.exp_mol_len = None
        self.mean_mol_len = None
        self.max_mol_len = None
        self.min_mol_len = None
        self.median_mol_len = None
##        self.mol_mol_len = None
        self.sum_2d_lengths = None
        self.nx_2d = None
        self.lx_2d = None
        self.exp_2d_len = None
        self.mean_2d_len = None
        self.median_2d_len = None
        self.min_2d_len = None
        self.max_2d_len = None
        self.sum_temp_lengths = None
        self.nx_temp = None
        self.lx_temp = None
        self.exp_temp_len = None
        self.mean_temp_len = None
        self.median_temp_len = None
        self.min_temp_len = None
        self.max_temp_len = None
        self.sum_comp_lengths = None
        self.nx_comp = None
        self.lx_comp = None
        self.exp_comp_len = None
        self.mean_comp_len = None
        self.median_comp_len = None
        self.min_comp_len = None
        self.max_comp_len = None
##        self.fxn_header = ["n_molecules","n_temp_only","n_comp","n_no_2d","n_comp_has_2d","n_comp_no_2d","n_2d","pct_temp_only","pct_has_comp","pct_has_2d","pct_no_2d","pct_no_2d_because_template_only","pct_no_2d_that_have_complement","pct_with_comp_that_has_2d","pct_with_comp_that_do_not_have_2d","sum_molecule_lengths","n_mol_gt","sum_mol_gt","pct_mol_gt","pct_data_from_mol_gt","mol_nx","mol_lx","expected_mol_len","mean_mol_len","median_mol_len","min_mol_len","max_mol_len","sum_2d_lengths","2d_nx","2d_lx","expected_2d_len","mean_2d_len","median_2d_len","min_2d_len","max_2d_len","sum_temp_lengths","temp_nx","temp_lx","expected_temp_len","mean_temp_len","median_temp_len","min_temp_len","max_temp_len","sum_comp_lengths","comp_nx","comp_lx","expected_comp_len","mean_comp_len","median_comp_len","min_comp_len","max_comp_len"]
##        self.allfxn = [self.get_n_molecules,self.get_n_temp_only,self.get_n_comp,self.get_n_no_2d,self.get_n_comp_has_2d,self.get_n_comp_no_2d,self.get_n_2d,self.get_pct_temp_only,self.get_pct_has_comp,self.get_pct_has_2d,self.get_pct_no_2d,self.get_pct_no_2d_because_template_only,self.get_pct_no_2d_that_have_complement,self.get_pct_with_comp_that_has_2d,self.get_pct_with_comp_that_do_not_have_2d,self.get_sum_molecule_lengths,self.get_n_mol_gt,self.get_sum_mol_gt,self.get_pct_mol_gt,self.get_pct_data_from_mol_gt,self.get_mol_nx,self.get_mol_lx,self.get_expected_mol_len,self.get_mean_mol_len,self.get_median_mol_len,self.get_min_mol_len,self.get_max_mol_len,self.get_sum_2d_lengths,self.get_2d_nx,self.get_2d_lx,self.get_expected_2d_len,self.get_mean_2d_len,self.get_median_2d_len,self.get_min_2d_len,self.get_max_2d_len,self.get_sum_temp_lengths,self.get_temp_nx,self.get_temp_lx,self.get_expected_temp_len,self.get_mean_temp_len,self.get_median_temp_len,self.get_min_temp_len,self.get_max_temp_len,self.get_sum_comp_lengths,self.get_comp_nx,self.get_comp_lx,self.get_expected_comp_len,self.get_mean_comp_len,self.get_median_comp_len,self.get_min_comp_len,self.get_max_comp_len]
        self.fxn_header = ["n_molecules","n_temp","n_temp_only","n_comp","n_no_2d","n_comp_has_2d","n_comp_no_2d","n_2d","pct_temp_only","pct_has_comp","pct_has_2d","pct_no_2d","pct_no_2d_because_template_only","pct_no_2d_that_have_complement","pct_with_comp_that_has_2d","pct_with_comp_that_do_not_have_2d","sum_molecule_lengths","n_mol_gt","sum_mol_gt","pct_mol_gt","pct_data_from_mol_gt","mol_nx","mol_lx","expected_mol_len","mean_mol_len","median_mol_len","min_mol_len","q_of_min_mol_len","max_mol_len","q_of_max_mol_len","mean_mol_q","median_mol_q","min_mol_q","len_of_min_mol_q","max_mol_q","len_of_max_mol_q","sum_2d_lengths","n_2d_gt","sum_2d_gt","pct_2d_gt","pct_2ddata_from_2d_gt","2d_nx","2d_lx","expected_2d_len","mean_2d_len","median_2d_len","min_2d_len","q_of_min_2d_len","max_2d_len","q_of_max_2d_len","mean_2d_q","median_2d_q","min_2d_q","len_of_min_2d_q","max_2d_q","len_of_max_2d_q","sum_temp_lengths","n_temp_gt","sum_temp_gt","pct_temp_gt","pct_tempdata_from_temp_gt","temp_nx","temp_lx","expected_temp_len","mean_temp_len","median_temp_len","min_temp_len","q_of_min_temp_len","max_temp_len","q_of_max_temp_len","mean_temp_q","median_temp_q","min_temp_q","len_of_min_temp_q","max_temp_q","len_of_max_temp_q","sum_comp_lengths","n_comp_gt","sum_comp_gt","pct_comp_gt","pct_compdata_from_comp_gt","comp_nx","comp_lx","expected_comp_len","mean_comp_len","median_comp_len","min_comp_len","q_of_min_comp_len","max_comp_len","q_of_max_comp_len","mean_comp_q","median_comp_q","min_comp_q","len_of_min_comp_q","max_comp_q","len_of_max_comp_q"]
        self.allfxn = [self.get_n_molecules,self.get_n_temp,self.get_n_temp_only,self.get_n_comp,self.get_n_no_2d,self.get_n_comp_has_2d,self.get_n_comp_no_2d,self.get_n_2d,self.get_pct_temp_only,self.get_pct_has_comp,self.get_pct_has_2d,self.get_pct_no_2d,self.get_pct_no_2d_because_template_only,self.get_pct_no_2d_that_have_complement,self.get_pct_with_comp_that_has_2d,self.get_pct_with_comp_that_do_not_have_2d,self.get_sum_molecule_lengths,self.get_n_mol_gt,self.get_sum_mol_gt,self.get_pct_mol_gt,self.get_pct_data_from_mol_gt,self.get_mol_nx,self.get_mol_lx,self.get_expected_mol_len,self.get_mean_mol_len,self.get_median_mol_len,self.get_min_mol_len,self.get_q_of_min_mol_len,self.get_max_mol_len,self.get_q_of_max_mol_len,self.get_mean_mol_q,self.get_median_mol_q,self.get_min_mol_q,self.get_len_of_min_mol_q,self.get_max_mol_q,self.get_len_of_max_mol_q,self.get_sum_2d_lengths,self.get_n_2d_gt,self.get_sum_2d_gt,self.get_pct_2d_gt,self.get_pct_2ddata_from_2d_gt,self.get_2d_nx,self.get_2d_lx,self.get_expected_2d_len,self.get_mean_2d_len,self.get_median_2d_len,self.get_min_2d_len,self.get_q_of_min_2d_len,self.get_max_2d_len,self.get_q_of_max_2d_len,self.get_mean_2d_q,self.get_median_2d_q,self.get_min_2d_q,self.get_len_of_min_2d_q,self.get_max_2d_q,self.get_len_of_max_2d_q,self.get_sum_temp_lengths,self.get_n_temp_gt,self.get_sum_temp_gt,self.get_pct_temp_gt,self.get_pct_tempdata_from_temp_gt,self.get_temp_nx,self.get_temp_lx,self.get_expected_temp_len,self.get_mean_temp_len,self.get_median_temp_len,self.get_min_temp_len,self.get_q_of_min_temp_len,self.get_max_temp_len,self.get_q_of_max_temp_len,self.get_mean_temp_q,self.get_median_temp_q,self.get_min_temp_q,self.get_len_of_min_temp_q,self.get_max_temp_q,self.get_len_of_max_temp_q,self.get_sum_comp_lengths,self.get_n_comp_gt,self.get_sum_comp_gt,self.get_pct_comp_gt,self.get_pct_compdata_from_comp_gt,self.get_comp_nx,self.get_comp_lx,self.get_expected_comp_len,self.get_mean_comp_len,self.get_median_comp_len,self.get_min_comp_len,self.get_q_of_min_comp_len,self.get_max_comp_len,self.get_q_of_max_comp_len,self.get_mean_comp_q,self.get_median_comp_q,self.get_min_comp_q,self.get_len_of_min_comp_q,self.get_max_comp_q,self.get_len_of_max_comp_q]
        self.bool_mol_len_gt = {}
        self.bool_mol_q_ge = {}
        self.n_mol_gt = {}
        self.sum_mol_gt = {}
        self.bool_2d_len_gt = {}
        self.bool_2d_q_ge = {}
        self.n_2d_gt = {}
        self.sum_2d_gt = {}
        self.bool_temp_len_gt = {}
        self.bool_temp_q_ge = {}
        self.n_temp_gt = {}
        self.sum_temp_gt = {}
        self.bool_comp_len_gt = {}
        self.bool_comp_q_ge = {}
        self.n_comp_gt = {}
        self.sum_comp_gt = {}
        #Q
        self.mean_mol_q = None
        self.median_mol_q = None
        self.min_mol_q = None
        self.max_mol_q = None
        self.mean_2d_q = None
        self.median_2d_q = None
        self.min_2d_q = None
        self.max_2d_q = None
        self.mean_temp_q = None
        self.median_temp_q = None
        self.min_temp_q = None
        self.max_temp_q = None
        self.mean_comp_q = None
        self.median_comp_q = None
        self.min_comp_q = None
        self.max_comp_q = None
        
    def get_n_molecules(self):
        if self.n_molecules is None:
            self.n_molecules = len(self.df['name'])
        return self.n_molecules

    def get_n_temp(self):
        ##At least for now, since all molecules with at least 1 read have a template read, numtemp = num.mol
        return self.get_n_molecules()
    
    def get_n_temp_only(self):
        if self.n_temp_only is None:
            self.n_temp_only = sum(self.temponly)
        return self.n_temp_only
    
    def get_n_comp(self):
        if self.n_comp is None:
            self.n_comp = sum(self.hascomp)
        return self.n_comp
    
    def get_n_no_2d(self):
        if self.n_no_2d is None:
            self.n_no_2d = sum(self.no2d)
        return self.n_no_2d
    
    def get_n_comp_has_2d(self): #should be same as number with 2d
        if self.n_comp_has_2d is None:
            self.n_comp_has_2d = sum(self.hascomp[self.has2d])
        return self.n_comp_has_2d
    
    def get_n_comp_no_2d(self):
        if self.n_comp_no_2d is None:
            self.n_comp_no_2d = sum(self.hascomp[self.no2d])
        return self.n_comp_no_2d
    
    def get_n_2d(self):
        if self.n_2d is None:
            self.n_2d = sum(self.has2d)
        return self.n_2d
    
    def get_pct_temp_only(self):
        return 100.0*self.get_n_temp_only()/self.get_n_molecules()

    def get_pct_has_comp(self):
        return 100.0*self.get_n_comp()/self.get_n_molecules()
    
    def get_pct_has_2d(self):
        return 100.0*self.get_n_2d()/self.get_n_molecules()
    
    def get_pct_no_2d(self):
        return 100.0*self.get_n_no_2d()/self.get_n_molecules()
    
    def get_pct_no_2d_because_template_only(self):
        return 100.0*self.get_n_temp_only()/self.get_n_no_2d()
    
    def get_pct_no_2d_that_have_complement(self):
        return 100.0*self.get_n_comp_no_2d()/self.get_n_no_2d()
    
    def get_pct_with_comp_that_has_2d(self):
        return 100.0*self.get_n_comp_has_2d()/self.get_n_comp()
    
    def get_pct_with_comp_that_do_not_have_2d(self):
        return 100.0*self.get_n_comp_no_2d()/self.get_n_comp()

############################################################
##'''Molecule'''
############################################################ 
    def get_sum_molecule_lengths(self):
        if self.sum_molecule_lengths is None:
            self.sum_molecule_lengths = sum(self.df['mol_len'])
        return self.sum_molecule_lengths

    def get_n_mol_gt(self, size=50e3, q=0):
        try:
            return self.n_mol_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='mol', force_redo = True)
            return self.n_mol_gt[size][q]

    def get_sum_mol_gt(self, size=50e3, q=0):
        try:
            return self.sum_mol_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='mol', force_redo = True)
            return self.sum_mol_gt[size][q]

    def get_pct_mol_gt(self, size=50e3,q=0):
        return 100.0*self.get_n_mol_gt(size,q)/self.get_n_molecules()

    def get_pct_data_from_mol_gt(self, size=50e3, q=0):
        return 100.0*self.get_sum_mol_gt(size,q)/self.get_sum_molecule_lengths()

    def get_mol_nx(self, x=[25,50,75], force_redo=False):
        if self.nx_mol is None or force_redo:
            self.nx_mol, self.lx_mol = self.NX(self.df['mol_len'], x=x, G=self.get_sum_molecule_lengths())
        return self.nx_mol

    def get_mol_lx(self, x=[25,50,75]):
        if self.nx_mol is None:
            self.get_mol_nx()
        return self.lx_mol

    def get_expected_mol_len(self):
        if self.exp_mol_len is None:
            self.exp_mol_len = self.e_size(self.df['mol_len'], G = self.get_sum_molecule_lengths())
        return self.exp_mol_len
        

    def get_mean_mol_len(self):
        if self.mean_mol_len is None:
            self.mean_mol_len = np.mean(self.df['mol_len'])
        return self.mean_mol_len

    def get_median_mol_len(self):
        if self.median_mol_len is None:
            self.median_mol_len = np.median(self.df['mol_len'])
        return self.median_mol_len


    def get_min_mol_len(self):
        if self.min_mol_len is None:
            self.min_mol_idx, self.min_mol_len = self.get_min(self.df['mol_len'])
            self.q_of_min_mol_len = self.df['mol_q'][self.min_mol_idx]
        return self.min_mol_len

    def get_q_of_min_mol_len(self):
        if self.min_mol_len is None:
            self.get_min_mol_len()
        return self.q_of_min_mol_len 
    
    def get_max_mol_len(self):
        if self.max_mol_len is None:
            self.max_mol_idx, self.max_mol_len = self.get_max(self.df['mol_len'])
            self.q_of_max_mol_len = self.df['mol_q'][self.max_mol_idx]
        return self.max_mol_len

    def get_q_of_max_mol_len(self):
        if self.max_mol_len is None:
            self.get_max_mol_len()
        return self.q_of_max_mol_len 

########quality
    def get_mean_mol_q(self):
        if self.mean_mol_q is None:
            self.mean_mol_q = np.mean(self.df['mol_q'])
        return self.mean_mol_q

    def get_median_mol_q(self):
        if self.median_mol_q is None:
            self.median_mol_q = np.median(self.df['mol_q'])
        return self.median_mol_q


    def get_min_mol_q(self):
        if self.min_mol_q is None:
            self.min_mol_q_idx, self.min_mol_q = self.get_min(self.df['mol_q'])
            self.len_of_min_mol_q = self.df['mol_len'][self.min_mol_q_idx]
        return self.min_mol_q

    def get_len_of_min_mol_q(self):
        if self.min_mol_q is None:
            self.get_min_mol_q()
        return self.len_of_min_mol_q
    
    def get_max_mol_q(self):
        if self.max_mol_q is None:
            self.max_mol_q_idx, self.max_mol_q = self.get_max(self.df['mol_q'])
            self.len_of_max_mol_q = self.df['mol_len'][self.max_mol_q_idx]
        return self.max_mol_q

    def get_len_of_max_mol_q(self):
        if self.max_mol_q is None:
            self.get_max_mol_q()
        return self.len_of_max_mol_q


############################################################
##'''2d'''
############################################################ 
    def get_sum_2d_lengths(self):
        if self.sum_2d_lengths is None:
##            #self.df['2d_len'].notnull()
            self.sum_2d_lengths = np.sum(self.df['2d_len'])
##            print sum(self.df['2d_len'][self.df['2d_len'].notnull()]) == np.sum(self.df['2d_len'][self.df['2d_len'].notnull()]), "SUM2d"
        return self.sum_2d_lengths

    def get_n_2d_gt(self, size=50e3, q=0):
        try:
            return self.n_2d_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='2d', force_redo = True)
            return self.n_2d_gt[size][q]

    def get_sum_2d_gt(self, size=50e3, q=0):
        try:
            return self.sum_2d_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='2d', force_redo = True)
            return self.sum_2d_gt[size][q]

    def get_pct_2d_gt(self, size=50e3,q=0):
        return 100.0*self.get_n_2d_gt(size,q)/self.get_n_2d()

    def get_pct_2ddata_from_2d_gt(self, size=50e3, q=0):
        return 100.0*self.get_sum_2d_gt(size,q)/self.get_sum_2d_lengths()

    def get_2d_nx(self, x=[25,50,75], force_redo=False):
        if self.nx_2d is None or force_redo:
            l = self.df['2d_len'][self.df['2d_len'].notnull()]
            self.nx_2d, self.lx_2d = self.NX(l, x=x, G=self.get_sum_2d_lengths())
        return self.nx_2d

    def get_2d_lx(self, x=[25,50,75]):
        if self.nx_2d is None:
            l = self.df['2d_len'][self.df['2d_len'].notnull()]
            self.nx_2d, self.lx_2d = self.NX(l, x=x, G=self.get_sum_2d_lengths())
        return self.lx_2d

    def get_expected_2d_len(self):
        if self.exp_2d_len is None:
            l = self.df['2d_len'][self.df['2d_len'].notnull()]
            self.exp_2d_len = self.e_size(l, G = self.get_sum_2d_lengths())
        return self.exp_2d_len
        

    def get_mean_2d_len(self):
        if self.mean_2d_len is None:
            self.mean_2d_len = np.mean(self.df['2d_len'])
        return self.mean_2d_len

    def get_median_2d_len(self):
        if self.median_2d_len is None:
            l = self.df['2d_len'][self.df['2d_len'].notnull()]
            self.median_2d_len = np.median(l)
        return self.median_2d_len


    def get_min_2d_len(self):
        if self.min_2d_len is None:
            self.min_2d_idx, self.min_2d_len = self.get_min(self.df['2d_len'])
            self.q_of_min_2d_len = self.df['2d_q'][self.min_2d_idx]
        return self.min_2d_len

    def get_q_of_min_2d_len(self):
        if self.min_2d_len is None:
            self.get_min_2d_len()
        return self.q_of_min_2d_len 
    
    def get_max_2d_len(self):
        if self.max_2d_len is None:
            self.max_2d_idx, self.max_2d_len = self.get_max(self.df['2d_len'])
            self.q_of_max_2d_len = self.df['2d_q'][self.max_2d_idx]
        return self.max_2d_len

    def get_q_of_max_2d_len(self):
        if self.max_2d_len is None:
            self.get_max_2d_len()
        return self.q_of_max_2d_len 

########quality
    def get_mean_2d_q(self):
        if self.mean_2d_q is None:
            self.mean_2d_q = np.mean(self.df['2d_q'])
        return self.mean_2d_q

    def get_median_2d_q(self):
        if self.median_2d_q is None:
            l = self.df['2d_q'][self.df['2d_q'].notnull()]
            self.median_2d_q = np.median(l)
        return self.median_2d_q


    def get_min_2d_q(self):
        if self.min_2d_q is None:
            self.min_2d_q_idx, self.min_2d_q = self.get_min(self.df['2d_q'])
            self.len_of_min_2d_q = self.df['2d_len'][self.min_2d_q_idx]
        return self.min_2d_q

    def get_len_of_min_2d_q(self):
        if self.min_2d_q is None:
            self.get_min_2d_q()
        return self.len_of_min_2d_q
    
    def get_max_2d_q(self):
        if self.max_2d_q is None:
            self.max_2d_q_idx, self.max_2d_q = self.get_max(self.df['2d_q'])
            self.len_of_max_2d_q = self.df['2d_len'][self.max_2d_q_idx]
        return self.max_2d_q

    def get_len_of_max_2d_q(self):
        if self.max_2d_q is None:
            self.get_max_2d_q()
        return self.len_of_max_2d_q

############################################################
##'''Temp'''
############################################################ 
    def get_sum_temp_lengths(self):
        if self.sum_temp_lengths is None:
            self.sum_temp_lengths = np.sum(self.df['temp_len'])
        return self.sum_temp_lengths

    def get_n_temp_gt(self, size=50e3, q=0):
        try:
            return self.n_temp_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='temp', force_redo = True)
            return self.n_temp_gt[size][q]

    def get_sum_temp_gt(self, size=50e3, q=0):
        try:
            return self.sum_temp_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='temp', force_redo = True)
            return self.sum_temp_gt[size][q]

    def get_pct_temp_gt(self, size=50e3,q=0):
        return 100.0*self.get_n_temp_gt(size,q)/self.get_n_temp()

    def get_pct_tempdata_from_temp_gt(self, size=50e3, q=0):
        return 100.0*self.get_sum_temp_gt(size,q)/self.get_sum_temp_lengths()   

    def get_temp_nx(self, x=[25,50,75], force_redo=False):
        if self.nx_temp is None or force_redo:
            l = self.df['temp_len'][self.df['temp_len'].notnull()]
            self.nx_temp, self.lx_temp = self.NX(l, x=x, G=self.get_sum_temp_lengths())
        return self.nx_temp

    def get_temp_lx(self,x=[25,50,75]):
        if self.nx_temp is None:
            l = self.df['temp_len'][self.df['temp_len'].notnull()]
            self.nx_temp, self.lx_temp = self.NX(l, x=x, G=self.get_sum_temp_lengths())
        return self.lx_temp

    def get_expected_temp_len(self):
        if self.exp_temp_len is None:
            l = self.df['temp_len'][self.df['temp_len'].notnull()]
            self.exp_temp_len = self.e_size(l, G = self.get_sum_temp_lengths())
        return self.exp_temp_len
        

    def get_mean_temp_len(self):
        if self.mean_temp_len is None:
            self.mean_temp_len = np.mean(self.df['temp_len'])
        return self.mean_temp_len

    def get_median_temp_len(self):
        if self.median_temp_len is None:
            l = self.df['temp_len'][self.df['temp_len'].notnull()]
            self.median_temp_len = np.median(l)
        return self.median_temp_len


    def get_min_temp_len(self):
        if self.min_temp_len is None:
            self.min_temp_idx, self.min_temp_len = self.get_min(self.df['temp_len'])
            self.q_of_min_temp_len = self.df['temp_q'][self.min_temp_idx]
        return self.min_temp_len
    
    def get_q_of_min_temp_len(self):
        if self.min_temp_len is None:
            self.get_min_temp_len()
        return self.q_of_min_temp_len 

    def get_max_temp_len(self):
        if self.max_temp_len is None:
            self.max_temp_idx, self.max_temp_len = self.get_max(self.df['temp_len'])
            self.q_of_max_temp_len = self.df['temp_q'][self.max_temp_idx]
        return self.max_temp_len

    def get_q_of_max_temp_len(self):
        if self.max_temp_len is None:
            self.get_max_temp_len()
        return self.q_of_max_temp_len

########quality
    def get_mean_temp_q(self):
        if self.mean_temp_q is None:
            self.mean_temp_q = np.mean(self.df['temp_q'])
        return self.mean_temp_q

    def get_median_temp_q(self):
        if self.median_temp_q is None:
            l = self.df['temp_q'][self.df['temp_q'].notnull()]
            self.median_temp_q = np.median(l)
        return self.median_temp_q


    def get_min_temp_q(self):
        if self.min_temp_q is None:
            self.min_temp_q_idx, self.min_temp_q = self.get_min(self.df['temp_q'])
            self.len_of_min_temp_q = self.df['temp_len'][self.min_temp_q_idx]
        return self.min_temp_q

    def get_len_of_min_temp_q(self):
        if self.min_temp_q is None:
            self.get_min_temp_q()
        return self.len_of_min_temp_q
    
    def get_max_temp_q(self):
        if self.max_temp_q is None:
            self.max_temp_q_idx, self.max_temp_q = self.get_max(self.df['temp_q'])
            self.len_of_max_temp_q = self.df['temp_len'][self.max_temp_q_idx]
        return self.max_temp_q

    def get_len_of_max_temp_q(self):
        if self.max_temp_q is None:
            self.get_max_temp_q()
        return self.len_of_max_temp_q
    
############################################################
##'''Comp'''
############################################################    
    def get_sum_comp_lengths(self):
        if self.sum_comp_lengths is None:
            self.sum_comp_lengths = np.sum(self.df['comp_len'])
        return self.sum_comp_lengths

    def get_n_comp_gt(self, size=50e3, q=0):
        try:
            return self.n_comp_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='comp', force_redo = True)
            return self.n_comp_gt[size][q]

    def get_sum_comp_gt(self, size=50e3, q=0):
        try:
            return self.sum_comp_gt[size][q]
        except:
            self._define_bools(size=size, q=q, read_type='comp', force_redo = True)
            return self.sum_comp_gt[size][q]

    def get_pct_comp_gt(self, size=50e3,q=0):
        return 100.0*self.get_n_comp_gt(size,q)/self.get_n_comp()

    def get_pct_compdata_from_comp_gt(self, size=50e3, q=0):
        return 100.0*self.get_sum_comp_gt(size,q)/self.get_sum_comp_lengths()

    def get_comp_nx(self, x=[25,50,75], force_redo=False):
        if self.nx_comp is None or force_redo:
            l = self.df['comp_len'][self.df['comp_len'].notnull()]
            self.nx_comp, self.lx_comp = self.NX(l, x=x, G=self.get_sum_comp_lengths())
        return self.nx_comp

    def get_comp_lx(self,x=[25,50,75]):
        if self.nx_comp is None:
            l = self.df['comp_len'][self.df['comp_len'].notnull()]
            self.nx_comp, self.lx_comp = self.NX(l, x=x, G=self.get_sum_comp_lengths())
        return self.lx_comp

    def get_expected_comp_len(self):
        if self.exp_comp_len is None:
            l = self.df['comp_len'][self.df['comp_len'].notnull()]
            self.exp_comp_len = self.e_size(l, G = self.get_sum_comp_lengths())
        return self.exp_comp_len
        

    def get_mean_comp_len(self):
        if self.mean_comp_len is None:
            self.mean_comp_len = np.mean(self.df['comp_len'])
        return self.mean_comp_len

    def get_median_comp_len(self):
        if self.median_comp_len is None:
            l = self.df['comp_len'][self.df['comp_len'].notnull()]
            self.median_comp_len = np.median(l)
        return self.median_comp_len


    def get_min_comp_len(self):
        if self.min_comp_len is None:
            self.min_comp_idx, self.min_comp_len = self.get_min(self.df['comp_len'])
            self.q_of_min_comp_len = self.df['comp_q'][self.min_comp_idx]
        return self.min_comp_len
    
    def get_q_of_min_comp_len(self):
        if self.min_comp_len is None:
            self.get_min_comp_len()
        return self.q_of_min_comp_len
    
    def get_max_comp_len(self):
        if self.max_comp_len is None:
            self.max_comp_idx, self.max_comp_len = self.get_max(self.df['comp_len'])
            self.q_of_max_comp_len = self.df['comp_q'][self.max_comp_idx]
        return self.max_comp_len
    
    def get_q_of_max_comp_len(self):
        if self.max_comp_len is None:
            self.get_max_comp_len()
        return self.q_of_max_comp_len 

########quality
    def get_mean_comp_q(self):
        if self.mean_comp_q is None:
            self.mean_comp_q = np.mean(self.df['comp_q'])
        return self.mean_comp_q

    def get_median_comp_q(self):
        if self.median_comp_q is None:
            l = self.df['comp_q'][self.df['comp_q'].notnull()]
            self.median_comp_q = np.median(l)
        return self.median_comp_q


    def get_min_comp_q(self):
        if self.min_comp_q is None:
            self.min_comp_q_idx, self.min_comp_q = self.get_min(self.df['comp_q'])
            self.len_of_min_comp_q = self.df['comp_len'][self.min_comp_q_idx]
        return self.min_comp_q

    def get_len_of_min_comp_q(self):
        if self.min_comp_q is None:
            self.get_min_comp_q()
        return self.len_of_min_comp_q
    
    def get_max_comp_q(self):
        if self.max_comp_q is None:
            self.max_comp_q_idx, self.max_comp_q = self.get_max(self.df['comp_q'])
            self.len_of_max_comp_q = self.df['comp_len'][self.max_comp_q_idx]
        return self.max_comp_q

    def get_len_of_max_comp_q(self):
        if self.max_comp_q is None:
            self.get_max_comp_q()
        return self.len_of_max_comp_q

############################################
##'''Top N Longest'''
############################################
    def N_longest_2d(self, N=10, q=0):
        self._define_bools(size=0, q=q, read_type='2d', force_redo = False, calc_gt=False)
        return self.df[self.bool_2d_q_ge[q]].sort(["2d_len","2d_q"], ascending=False)[:N]

    def N_longest_temp(self, N=10, q=0):
        self._define_bools(size=0, q=q, read_type='temp', force_redo = False, calc_gt=False)
        return self.df[self.bool_temp_q_ge[q]].sort(["temp_len","temp_q"], ascending=False)[:N]

    def N_longest_comp(self, N=10, q=0):
        self._define_bools(size=0, q=q, read_type='comp', force_redo = False, calc_gt=False)
        return self.df[self.bool_comp_q_ge[q]].sort(["comp_len","comp_q"], ascending=False)[:N]

    def N_longest_dict(self, read_type="2d", N=10, q=0, columns=["temp_len","temp_q","name"]):
        nl = {"2d":self.N_longest_2d, "temp":self.N_longest_temp, "comp":self.N_longest_comp}
        df = nl[read_type](N,q)[columns]
        out = {}
        for col in columns:
            out[col] = list(df[col])
        return out

    def print_N_longest(self, read_type="2d", N=10, q=0, rank=True, columns=None, add_read_type=True):
        if columns is None:
            columns = [read_type+"_len",read_type+"_q","name"]
        rt=""
        if add_read_type:
            rt = read_type+"\t"
        out = self.N_longest_dict(read_type, N, q, columns)
        if rank:
            for i in range(len(out[columns[0]])):
                print rt+("\t").join([str(e) for e in [i+1] + [out[j][i] for j in columns]])
        else:
            for i in range(len(out[columns[0]])):
                print rt+("\t").join([str(e) for e in [out[j][i] for j in columns]])


############################################
##'''tools'''
############################################
    def get_min(self,l):
        idx = l.idxmin()
        m = l[idx]
        return idx, m

    def get_max(self,l):
        idx = l.idxmax()
        m = l[idx]
        return idx, m

    def get_max_len_after_filter_q(self, q=0, read_type="2d"):
        self._define_bools(size=0, q=q, read_type=read_type, force_redo = False, calc_gt=False)
        Q = {"mol":self.bool_mol_q_ge, "2d":self.bool_2d_q_ge, "temp":self.bool_temp_q_ge, "comp":self.bool_comp_q_ge}
        r_len = read_type + "_len"
        r_q = read_type + "_q"
        idx, max_len = self.get_max(self.df[r_len][Q[read_type][q]])
        max_len_q = self.df[r_q][Q[read_type][q]][idx]
        return max_len, max_len_q

    def get_max_q(self, read_type="mol"):
        Q = {"mol":self.get_max_mol_q, "2d":self.get_max_2d_q, "temp":self.get_max_temp_q, "comp":self.get_max_comp_q}
        return Q[read_type]()
        

    def NX(self, l, x=[25,50,75], G=False):
        """
        Returns NX for all x for a list of numbers l.
        Default: N25, N50, N75
        Assumes all values in list x are between 0 and 100.
        Interpretation: When NX = NX_value, X% of data (in bp) is contained in reads at least NX_value bp long.
        """
	if G:
            l = sorted(list(l[:]))
            x = sorted(list(x[:]))
            nxsum = 0
            L = 0
            nxvalues = {e:0 for e in x}
            lxvalues = {e:0 for e in x}
            for e in x:
                    xpct = G*e/100.0
                    while nxsum < xpct and l:
                            nxsum += l[-1]
                            L += 1
                            lastsize = l.pop()
                    nxvalues[e] = lastsize
                    lxvalues[e] = L
            return nxvalues, lxvalues

    def e_size(self, l,G=False):
        if not G:
            G = sum(l)
        return sum([e**2 for e in list(l)])/float(G)
        
    def n_lengths_gt(self, l, size):
        #l from df
        ans = l > size
        return sum(ans)

    def pct_sum_from_lengths_gt(self, l, size):
        #l from df
        ans = l > size
        return sum(l[ans])

    def summarize(self, fxn_header, allfxn, x=[25,50,75], gt=[10e3, 50e3, 100e3], qmol=[4.5,6,7,8,10], q2d=[9,10,11,12], q1d=[3.5,4.5,6]):
        ##TODO:
        ## Add in longest with Q > x,y,z; for 2d Q(9,10,13); 1d Q(3.5,4.5,6)
        ## Add in n mol/2d/temp/comp gt L and gt Q --- though for mol do L(10)Q(3.5,4.5,6,7,8,10); 2d do L(10)Q(9,10,11,12); 1d do L(10)Q(3.5,4.5,6)
        summary = []
        summary_header = []
        x = sorted(x)
        for i in range(len(allfxn)):
            f = allfxn[i]
            h = fxn_header[i]
            if h[-2:] in ["nx", "lx"]:
                x = f()
                for e in x:
                    hh = h[:-1] + str(e)
                    xx = x[e]
                    summary_header.append(hh)
                    summary.append(xx)
##                    print ("\t").join([hh,str(xx)])
            elif h[-2:] == "gt":
                for e in gt:
                    hh = h+"_"+str(int(e/1000))+"kb"
                    xx = f(size = e)
                    summary_header.append(hh)
                    summary.append(xx)
##                    print ("\t").join([hh,str(xx)])
                if "mol" in h:
                    for e in qmol:
                        hh = h+"_"+str(int(gt[0]/1000))+"kb_ge_Q"+str(e)
                        xx = f(size = gt[0], q=e)
                        summary_header.append(hh)
                        summary.append(xx)
                elif "2d" in h:
                    for e in q2d:
                        hh = h+"_"+str(int(gt[0]/1000))+"kb_ge_Q"+str(e)
                        xx = f(size = gt[0], q=e)
                        summary_header.append(hh)
                        summary.append(xx)
                elif "temp" in h or "comp" in h:
                    for e in q1d:
                        hh = h+"_"+str(int(gt[0]/1000))+"kb_ge_Q"+str(e)
                        xx = f(size = gt[0], q=e)
                        summary_header.append(hh)
                        summary.append(xx)
            else:
                x = f()
                summary_header.append(h)
                summary.append(x)
##                print ("\t").join([h,str(x)])
                if "max" in h and "len" in h and "q_of" not in h:
                    if "mol" in h:
                        for e in qmol:
                            if e <= self.get_max_q("mol"):
                                max_len, max_len_q = self.get_max_len_after_filter_q(q=e, read_type="mol")
                            else:
                                max_len, max_len_q = "-", "-"
                            hh = h + "_Q_ge_"+str(e)
                            summary_header.append(hh)
                            summary.append(max_len)
                            hhh = "q_of_" + hh
                            summary_header.append(hhh)
                            summary.append(max_len_q)
                            
                    elif "2d" in h:
                        for e in q2d:
                            if e <= self.get_max_q("2d"):
                                max_len, max_len_q = self.get_max_len_after_filter_q(q=e, read_type="2d")
                            else:
                                max_len, max_len_q = "-", "-"
                            hh = h + "_Q_ge_"+str(e)
                            summary_header.append(hh)
                            summary.append(max_len)
                            hhh = "q_of_" + hh
                            summary_header.append(hhh)
                            summary.append(max_len_q)
                    if "temp" in h or "comp" in h:
                        if "temp" in h:
                            rt="temp"
                        else:
                            rt="comp"
                        for e in q1d:
                            if e <= self.get_max_q(rt):
                                max_len, max_len_q = self.get_max_len_after_filter_q(q=e, read_type=rt)
                            else:
                                max_len, max_len_q = "-", "-"
                            hh = h + "_Q_ge_"+str(e)
                            summary_header.append(hh)
                            summary.append(max_len)
                            hhh = "q_of_" + hh
                            summary_header.append(hhh)
                            summary.append(max_len_q)
                        
        return summary_header, summary

    def standard_summary(self, x=[25,50,75], gt=[10e3, 50e3, 100e3]):
        return self.summarize(fxn_header=self.fxn_header, allfxn=self.allfxn, x=x, gt=gt)


    def read_type_summary(self):
        return self.summarize(fxn_header=self.fxn_header[:15], allfxn=self.allfxn[:15])
 
    def _define_bools(self, size=0, q=0, read_type='mol', force_redo=False, calc_gt=True):
        except_len = False
        except_q = False
        if read_type == 'mol':
            try:
                self.bool_mol_len_gt[size]
            except:
                except_len = True
                self.bool_mol_len_gt[size] = self.df['mol_len'] > size
            try:
                self.bool_mol_q_ge[q]
            except:
                except_q = True
                self.bool_mol_q_ge[q] = self.df['mol_q'] >= q
            if calc_gt and (except_len or except_q or force_redo):
                try:
                    self.n_mol_gt[size]
                except:
                    self.n_mol_gt[size] = {}
                try:
                    self.sum_mol_gt[size]
                except:
                    self.sum_mol_gt[size] = {}
                ans = self.df['mol_len'][self.bool_mol_len_gt[size]][self.bool_mol_q_ge[q]]
                self.n_mol_gt[size][q] = len(ans)
                self.sum_mol_gt[size][q] = sum(ans)
                
        elif read_type == '2d':
            try:
                self.bool_2d_len_gt[size]
            except:
                except_len = True
                self.bool_2d_len_gt[size] = self.df['2d_len'] > size
            try:
                self.bool_2d_q_ge[q]
            except:
                except_q = True
                self.bool_2d_q_ge[q] = self.df['2d_q'] >= q
            if calc_gt and (except_len or except_q or force_redo):
                try:
                    self.n_2d_gt[size]
                except:
                    self.n_2d_gt[size] = {}
                try:
                    self.sum_2d_gt[size]
                except:
                    self.sum_2d_gt[size] = {}
                ans = self.df['2d_len'][self.bool_2d_len_gt[size]][self.bool_2d_q_ge[q]]
                self.n_2d_gt[size][q] = len(ans)
                self.sum_2d_gt[size][q] = sum(ans)
        elif read_type == 'temp':
            try:
                self.bool_temp_len_gt[size]
            except:
                except_len = True
                self.bool_temp_len_gt[size] = self.df['temp_len'] > size
            try:
                self.bool_temp_q_ge[q]
            except:
                except_q = True
                self.bool_temp_q_ge[q] = self.df['temp_q'] >= q
            if calc_gt and (except_len or except_q or force_redo):
                try:
                    self.n_temp_gt[size]
                except:
                    self.n_temp_gt[size] = {}
                try:
                    self.sum_temp_gt[size]
                except:
                    self.sum_temp_gt[size] = {}
                ans = self.df['temp_len'][self.bool_temp_len_gt[size]][self.bool_temp_q_ge[q]]
                self.n_temp_gt[size][q] = len(ans)
                self.sum_temp_gt[size][q] = sum(ans)
        elif read_type == 'comp':
            try:
                self.bool_comp_len_gt[size]
            except:
                except_len = True
                self.bool_comp_len_gt[size] = self.df['comp_len'] > size
            try:
                self.bool_comp_q_ge[q]
            except:
                except_q = True
                self.bool_comp_q_ge[q] = self.df['comp_q'] >= q
            if except_len or except_q or force_redo:
                try:
                    self.n_comp_gt[size]
                except:
                    self.n_comp_gt[size] = {}
                try:
                    self.sum_comp_gt[size]
                except:
                    self.sum_comp_gt[size] = {}
                ans = self.df['comp_len'][self.bool_comp_len_gt[size]][self.bool_comp_q_ge[q]]
                self.n_comp_gt[size][q] = len(ans)
                self.sum_comp_gt[size][q] = sum(ans)
