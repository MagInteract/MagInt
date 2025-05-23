from __future__ import print_function, absolute_import, division
from . import _hubbard_I
import f90wrap.runtime
import logging
import numpy

class Hubbard_I_Data(f90wrap.runtime.FortranModule):
    """
    Module hubbard_i_data
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 17-32
    
    """
    @f90wrap.runtime.register_class("hubbard_I.occup")
    class occup(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=occup)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 20-31
        
        """
        def __init__(self, handle=None):
            """
            self = Occup()
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                lines 20-31
            
            
            Returns
            -------
            this : Occup
            	Object to be constructed
            
            
            Automatically generated constructor for occup
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _hubbard_I.f90wrap_occup_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Occup
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                lines 20-31
            
            Parameters
            ----------
            this : Occup
            	Object to be destructed
            
            
            Automatically generated destructor for occup
            """
            if self._alloc:
                _hubbard_I.f90wrap_occup_finalise(this=self._handle)
        
        @property
        def ifdiag(self):
            """
            Element ifdiag ftype=logical pytype=bool
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 21
            
            """
            return _hubbard_I.f90wrap_occup__get__ifdiag(self._handle)
        
        @ifdiag.setter
        def ifdiag(self, ifdiag):
            _hubbard_I.f90wrap_occup__set__ifdiag(self._handle, ifdiag)
        
        @property
        def run(self):
            """
            Element run ftype=logical pytype=bool
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 21
            
            """
            return _hubbard_I.f90wrap_occup__get__run(self._handle)
        
        @run.setter
        def run(self, run):
            _hubbard_I.f90wrap_occup__set__run(self._handle, run)
        
        @property
        def ndeg(self):
            """
            Element ndeg ftype=integer  pytype=int
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 22
            
            """
            return _hubbard_I.f90wrap_occup__get__ndeg(self._handle)
        
        @ndeg.setter
        def ndeg(self, ndeg):
            _hubbard_I.f90wrap_occup__set__ndeg(self._handle, ndeg)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 23
            
            """
            return _hubbard_I.f90wrap_occup__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _hubbard_I.f90wrap_occup__set__n(self._handle, n)
        
        @property
        def hn(self):
            """
            Element hn ftype=complex(8) pytype=complex
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _hubbard_I.f90wrap_occup__array__hn(self._handle)
            if array_handle in self._arrays:
                hn = self._arrays[array_handle]
            else:
                hn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _hubbard_I.f90wrap_occup__array__hn)
                self._arrays[array_handle] = hn
            return hn
        
        @hn.setter
        def hn(self, hn):
            self.hn[...] = hn
        
        @property
        def en(self):
            """
            Element en ftype=real(8) pytype=float
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _hubbard_I.f90wrap_occup__array__en(self._handle)
            if array_handle in self._arrays:
                en = self._arrays[array_handle]
            else:
                en = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _hubbard_I.f90wrap_occup__array__en)
                self._arrays[array_handle] = en
            return en
        
        @en.setter
        def en(self, en):
            self.en[...] = en
        
        @property
        def st_n(self):
            """
            Element st_n ftype=integer pytype=int
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _hubbard_I.f90wrap_occup__array__st_n(self._handle)
            if array_handle in self._arrays:
                st_n = self._arrays[array_handle]
            else:
                st_n = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _hubbard_I.f90wrap_occup__array__st_n)
                self._arrays[array_handle] = st_n
            return st_n
        
        @st_n.setter
        def st_n(self, st_n):
            self.st_n[...] = st_n
        
        @property
        def narr(self):
            """
            Element narr ftype=integer pytype=int
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                line 31
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _hubbard_I.f90wrap_occup__array__narr(self._handle)
            if array_handle in self._arrays:
                narr = self._arrays[array_handle]
            else:
                narr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _hubbard_I.f90wrap_occup__array__narr)
                self._arrays[array_handle] = narr
            return narr
        
        @narr.setter
        def narr(self, narr):
            self.narr[...] = narr
        
        def __str__(self):
            ret = ['<occup>{\n']
            ret.append('    ifdiag : ')
            ret.append(repr(self.ifdiag))
            ret.append(',\n    run : ')
            ret.append(repr(self.run))
            ret.append(',\n    ndeg : ')
            ret.append(repr(self.ndeg))
            ret.append(',\n    n : ')
            ret.append(repr(self.n))
            ret.append(',\n    hn : ')
            ret.append(repr(self.hn))
            ret.append(',\n    en : ')
            ret.append(repr(self.en))
            ret.append(',\n    st_n : ')
            ret.append(repr(self.st_n))
            ret.append(',\n    narr : ')
            ret.append(repr(self.narr))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("hubbard_I.Occup_X100_Array")
    class Occup_X100_Array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=occup_x100_array)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 20-31
        
        super-type
        Automatically generated to handle derived type arrays as a new derived type
        """
        def __init__(self, handle=None):
            """
            self = Occup_X100_Array()
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                lines 20-31
            
            
            Returns
            -------
            this : Occup_X100_Array
            	Object to be constructed
            
            
            Automatically generated constructor for occup_x100_array
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _hubbard_I.f90wrap_occup_x100_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Occup_X100_Array
            
            
            Defined at \
                /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
                lines 20-31
            
            Parameters
            ----------
            this : Occup_X100_Array
            	Object to be destructed
            
            
            Automatically generated destructor for occup_x100_array
            """
            if self._alloc:
                _hubbard_I.f90wrap_occup_x100_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _hubbard_I.f90wrap_occup_x100_array__array_getitem__items,
                                            _hubbard_I.f90wrap_occup_x100_array__array_setitem__items,
                                            _hubbard_I.f90wrap_occup_x100_array__array_len__items,
                                            """
            Element items ftype=type(occup) pytype=Occup
            
            
            Defined at  line 0
            
            """, Hubbard_I_Data.occup)
            return self.items
        
        _dt_array_initialisers = [init_array_items]
        
    
    @property
    def b(self):
        """
        Element b ftype=real(8) pytype=float
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 19
        
        """
        return _hubbard_I.f90wrap_hubbard_i_data__get__b()
    
    @b.setter
    def b(self, b):
        _hubbard_I.f90wrap_hubbard_i_data__set__b(b)
    
    def __str__(self):
        ret = ['<hubbard_i_data>{\n']
        ret.append('    b : ')
        ret.append(repr(self.b))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

hubbard_i_data = Hubbard_I_Data()

class States_Gn(f90wrap.runtime.FortranModule):
    """
    Module states_gn
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 34-213
    
    """
    @staticmethod
    def c_dag_c_dag_cc(vec, m, m1, m2, m3):
        """
        c_dag_c_dag_cc = c_dag_c_dag_cc(vec, m, m1, m2, m3)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 45-85
        
        Parameters
        ----------
        vec : complex array
        m : int
        m1 : int
        m2 : int
        m3 : int
        
        Returns
        -------
        c_dag_c_dag_cc : complex array
        
        """
        c_dag_c_dag_cc = _hubbard_I.f90wrap_c_dag_c_dag_cc(vec=vec, m=m, m1=m1, m2=m2, \
            m3=m3)
        return c_dag_c_dag_cc
    
    @staticmethod
    def c_dag_c(vec, m, m1):
        """
        c_dag_c = c_dag_c(vec, m, m1)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 87-121
        
        Parameters
        ----------
        vec : complex array
        m : int
        m1 : int
        
        Returns
        -------
        c_dag_c : complex array
        
        """
        c_dag_c = _hubbard_I.f90wrap_c_dag_c(vec=vec, m=m, m1=m1)
        return c_dag_c
    
    @staticmethod
    def ovl(vec, vec1):
        """
        ovl = ovl(vec, vec1)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 170-182
        
        Parameters
        ----------
        vec : complex array
        vec1 : complex array
        
        Returns
        -------
        ovl : complex
        
        """
        ovl = _hubbard_I.f90wrap_ovl(vec=vec, vec1=vec1)
        return ovl
    
    @staticmethod
    def setup_aux_arrays(nat, nso):
        """
        n = setup_aux_arrays(nat, nso)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 184-206
        
        Parameters
        ----------
        nat : int
        nso : int
        
        Returns
        -------
        n : int
        
        """
        n = _hubbard_I.f90wrap_setup_aux_arrays(nat=nat, nso=nso)
        return n
    
    @staticmethod
    def free_aux_arrays():
        """
        free_aux_arrays()
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 208-212
        
        
        """
        _hubbard_I.f90wrap_free_aux_arrays()
    
    @staticmethod
    def _op_1el(vec, mat, nso, n):
        """
        op_1el = _op_1el(vec, mat, nso, n)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 151-168
        
        Parameters
        ----------
        vec : complex array
        mat : complex array
        nso : int
        n : int
        
        Returns
        -------
        op_1el : complex array
        
        """
        op_1el = _hubbard_I.f90wrap_op_1el(vec=vec, mat=mat, nso=nso, n=n)
        return op_1el
    
    @staticmethod
    def _op_2el(vec, mat, nso, n):
        """
        op_2el = _op_2el(vec, mat, nso, n)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 123-148
        
        Parameters
        ----------
        vec : complex array
        mat : complex array
        nso : int
        n : int
        
        Returns
        -------
        op_2el : complex array
        
        """
        op_2el = _hubbard_I.f90wrap_op_2el(vec=vec, mat=mat, nso=nso, n=n)
        return op_2el
    
    @staticmethod
    def op(*args, **kwargs):
        """
        op(*args, **kwargs)
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            lines 40-41
        
        Overloaded interface containing the following procedures:
          _op_1el
          _op_2el
        
        """
        for proc in [States_Gn._op_1el, States_Gn._op_2el]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
        
    
    @property
    def tol_mel(self):
        """
        Element tol_mel ftype=real(kind=8) pytype=float
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 35
        
        """
        return _hubbard_I.f90wrap_states_gn__get__tol_mel()
    
    @property
    def st_gn(self):
        """
        Element st_gn ftype=integer pytype=int
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hubbard_I.f90wrap_states_gn__array__st_gn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            st_gn = self._arrays[array_handle]
        else:
            st_gn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hubbard_I.f90wrap_states_gn__array__st_gn)
            self._arrays[array_handle] = st_gn
        return st_gn
    
    @st_gn.setter
    def st_gn(self, st_gn):
        self.st_gn[...] = st_gn
    
    @property
    def narr_g(self):
        """
        Element narr_g ftype=integer pytype=int
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hubbard_I.f90wrap_states_gn__array__narr_g(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            narr_g = self._arrays[array_handle]
        else:
            narr_g = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hubbard_I.f90wrap_states_gn__array__narr_g)
            self._arrays[array_handle] = narr_g
        return narr_g
    
    @narr_g.setter
    def narr_g(self, narr_g):
        self.narr_g[...] = narr_g
    
    @property
    def occ(self):
        """
        Element occ ftype=real(8) pytype=float
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hubbard_I.f90wrap_states_gn__array__occ(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            occ = self._arrays[array_handle]
        else:
            occ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hubbard_I.f90wrap_states_gn__array__occ)
            self._arrays[array_handle] = occ
        return occ
    
    @occ.setter
    def occ(self, occ):
        self.occ[...] = occ
    
    @property
    def arr(self):
        """
        Element arr ftype=integer pytype=int
        
        
        Defined at \
            /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
            line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _hubbard_I.f90wrap_states_gn__array__arr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            arr = self._arrays[array_handle]
        else:
            arr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _hubbard_I.f90wrap_states_gn__array__arr)
            self._arrays[array_handle] = arr
        return arr
    
    @arr.setter
    def arr(self, arr):
        self.arr[...] = arr
    
    def __str__(self):
        ret = ['<states_gn>{\n']
        ret.append('    tol_mel : ')
        ret.append(repr(self.tol_mel))
        ret.append(',\n    st_gn : ')
        ret.append(repr(self.st_gn))
        ret.append(',\n    narr_g : ')
        ret.append(repr(self.narr_g))
        ret.append(',\n    occ : ')
        ret.append(repr(self.occ))
        ret.append(',\n    arr : ')
        ret.append(repr(self.arr))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

states_gn = States_Gn()

def op_exp_val(op_mat, nat, expvals, nvec, nso):
    """
    op_exp_val(op_mat, nat, expvals, nvec, nso)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 215-238
    
    Parameters
    ----------
    op_mat : complex array
    nat : int
    expvals : complex array
    nvec : int
    nso : int
    
    """
    _hubbard_I.f90wrap_op_exp_val(op_mat=op_mat, nat=nat, expvals=expvals, \
        nvec=nvec, nso=nso)

def op_exp_val_invec(resvec, op_mat, invec, nat, nso, n):
    """
    expval = op_exp_val_invec(resvec, op_mat, invec, nat, nso, n)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 240-267
    
    Parameters
    ----------
    resvec : complex array
    op_mat : complex array
    invec : complex array
    nat : int
    nso : int
    n : int
    
    Returns
    -------
    expval : complex
    
    """
    expval = _hubbard_I.f90wrap_op_exp_val_invec(resvec=resvec, op_mat=op_mat, \
        invec=invec, nat=nat, nso=nso, n=n)
    return expval

def prn_vectors(nat, nso, n, st_num, tresh):
    """
    prn_vectors(nat, nso, n, st_num, tresh)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 269-301
    
    Parameters
    ----------
    nat : int
    nso : int
    n : int
    st_num : int
    tresh : float
    
    """
    _hubbard_I.f90wrap_prn_vectors(nat=nat, nso=nso, n=n, st_num=st_num, \
        tresh=tresh)

def prn_input_vector(nat, nso, n, vec, tresh):
    """
    prn_input_vector(nat, nso, n, vec, tresh)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 303-329
    
    Parameters
    ----------
    nat : int
    nso : int
    n : int
    vec : complex array
    tresh : float
    
    """
    _hubbard_I.f90wrap_prn_input_vector(nat=nat, nso=nso, n=n, vec=vec, tresh=tresh)

def get_vec(n, iv, vec):
    """
    get_vec(n, iv, vec)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 331-353
    
    Parameters
    ----------
    n : int
    iv : int
    vec : complex array
    
    """
    _hubbard_I.f90wrap_get_vec(n=n, iv=iv, vec=vec)

def mat_el_vecs(nat, nso, n, st_l, st_r, op_mat):
    """
    val = mat_el_vecs(nat, nso, n, st_l, st_r, op_mat)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 355-380
    
    Parameters
    ----------
    nat : int
    nso : int
    n : int
    st_l : complex array
    st_r : complex array
    op_mat : complex array
    
    Returns
    -------
    val : complex
    
    """
    val = _hubbard_I.f90wrap_mat_el_vecs(nat=nat, nso=nso, n=n, st_l=st_l, \
        st_r=st_r, op_mat=op_mat)
    return val

def mat_el_vecs_2el(nat, nso, n, st_l, st_r, op_mat):
    """
    val = mat_el_vecs_2el(nat, nso, n, st_l, st_r, op_mat)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 382-403
    
    Parameters
    ----------
    nat : int
    nso : int
    n : int
    st_l : complex array
    st_r : complex array
    op_mat : complex array
    
    Returns
    -------
    val : complex
    
    """
    val = _hubbard_I.f90wrap_mat_el_vecs_2el(nat=nat, nso=nso, n=n, st_l=st_l, \
        st_r=st_r, op_mat=op_mat)
    return val

def op_on_state(nat, nso, n, st_num, st_res, op_mat):
    """
    norm = op_on_state(nat, nso, n, st_num, st_res, op_mat)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 405-437
    
    Parameters
    ----------
    nat : int
    nso : int
    n : int
    st_num : int
    st_res : complex array
    op_mat : complex array
    
    Returns
    -------
    norm : float
    
    """
    norm = _hubbard_I.f90wrap_op_on_state(nat=nat, nso=nso, n=n, st_num=st_num, \
        st_res=st_res, op_mat=op_mat)
    return norm

def ops_on_state(st_res, st_init, op_mat, nat, nso, n, ntimes):
    """
    ops_on_state(st_res, st_init, op_mat, nat, nso, n, ntimes)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 439-462
    
    Parameters
    ----------
    st_res : complex array
    st_init : complex array
    op_mat : complex array
    nat : int
    nso : int
    n : int
    ntimes : int
    
    """
    _hubbard_I.f90wrap_ops_on_state(st_res=st_res, st_init=st_init, op_mat=op_mat, \
        nat=nat, nso=nso, n=n, ntimes=ntimes)

def fock_space_rmat(fs_rmat, rm, nat, nso, n):
    """
    fock_space_rmat(fs_rmat, rm, nat, nso, n)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 464-543
    
    Parameters
    ----------
    fs_rmat : complex array
    rm : complex array
    nat : int
    nso : int
    n : int
    
    """
    _hubbard_I.f90wrap_fock_space_rmat(fs_rmat=fs_rmat, rm=rm, nat=nat, nso=nso, \
        n=n)

def c_dag_to_set(vecs_in, coff_in, vecs_out, coff_out, n_in, ind_in, m, nso, \
    max_st):
    """
    n_out = c_dag_to_set(vecs_in, coff_in, vecs_out, coff_out, n_in, ind_in, m, nso, \
        max_st)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 545-581
    
    Parameters
    ----------
    vecs_in : int array
    coff_in : complex array
    vecs_out : int array
    coff_out : complex array
    n_in : int
    ind_in : int
    m : int
    nso : int
    max_st : int
    
    Returns
    -------
    n_out : int
    
    """
    n_out = _hubbard_I.f90wrap_c_dag_to_set(vecs_in=vecs_in, coff_in=coff_in, \
        vecs_out=vecs_out, coff_out=coff_out, n_in=n_in, ind_in=ind_in, m=m, \
        nso=nso, max_st=max_st)
    return n_out

def gf_hi_fullu_int(e0f, u, ummss, zmsb, nlm, iwmax, nmom, ns, temp, verbosity, \
    n_lev, calc_off_diag, remove_cf, ladbs, ladop, calcovl, nbas, mnat, stbas, \
    gf0, tail0, gf, tail, ovlmat):
    """
    atocc, atmag, z = gf_hi_fullu_int(e0f, u, ummss, zmsb, nlm, iwmax, nmom, ns, \
        temp, verbosity, n_lev, calc_off_diag, remove_cf, ladbs, ladop, calcovl, \
        nbas, mnat, stbas, gf0, tail0, gf, tail, ovlmat)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 583-956
    
    Parameters
    ----------
    e0f : complex array
    u : float array
    ummss : float array
    zmsb : complex array
    nlm : int
    iwmax : int
    nmom : int
    ns : int
    temp : float
    verbosity : int
    n_lev : int
    calc_off_diag : bool
    remove_cf : bool
    ladbs : bool
    ladop : complex array
    calcovl : bool
    nbas : int
    mnat : int
    stbas : complex array
    gf0 : complex array
    tail0 : complex array
    gf : complex array
    tail : complex array
    ovlmat : complex array
    
    Returns
    -------
    atocc : float
    atmag : float
    z : float
    
    """
    atocc, atmag, z = _hubbard_I.f90wrap_gf_hi_fullu_int(e0f=e0f, u=u, ummss=ummss, \
        zmsb=zmsb, nlm=nlm, iwmax=iwmax, nmom=nmom, ns=ns, temp=temp, \
        verbosity=verbosity, n_lev=n_lev, calc_off_diag=calc_off_diag, \
        remove_cf=remove_cf, ladbs=ladbs, ladop=ladop, calcovl=calcovl, nbas=nbas, \
        mnat=mnat, stbas=stbas, gf0=gf0, tail0=tail0, gf=gf, tail=tail, \
        ovlmat=ovlmat)
    return atocc, atmag, z

def factor(n):
    """
    factor = factor(n)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 958-968
    
    Parameters
    ----------
    n : int
    
    Returns
    -------
    factor : float
    
    """
    factor = _hubbard_I.f90wrap_factor(n=n)
    return factor

def diagh(h, e, e0f, u, st, arr, narr, nso, nstate, nat, n, verbosity):
    """
    diagh(h, e, e0f, u, st, arr, narr, nso, nstate, nat, n, verbosity)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 970-1061
    
    Parameters
    ----------
    h : complex array
    e : float array
    e0f : complex array
    u : float array
    st : int array
    arr : int array
    narr : int array
    nso : int
    nstate : int
    nat : int
    n : int
    verbosity : int
    
    """
    _hubbard_I.f90wrap_diagh(h=h, e=e, e0f=e0f, u=u, st=st, arr=arr, narr=narr, \
        nso=nso, nstate=nstate, nat=nat, n=n, verbosity=verbosity)

def add_to_gf_n(gf, tail, arr, nso, nmom, nat, iwmax, num, n_occ, zmsb, z, \
    eground, temp, left, right):
    """
    add_to_gf_n(gf, tail, arr, nso, nmom, nat, iwmax, num, n_occ, zmsb, z, eground, \
        temp, left, right)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 1063-1186
    
    Parameters
    ----------
    gf : complex array
    tail : complex array
    arr : int array
    nso : int
    nmom : int
    nat : int
    iwmax : int
    num : int
    n_occ : Occup_X100_Array
    	super-type
    
    zmsb : complex array
    z : float
    eground : float
    temp : float
    left : bool
    right : bool
    
    """
    _hubbard_I.f90wrap_add_to_gf_n(gf=gf, tail=tail, arr=arr, nso=nso, nmom=nmom, \
        nat=nat, iwmax=iwmax, num=num, n_occ=n_occ._handle, zmsb=zmsb, z=z, \
        eground=eground, temp=temp, left=left, right=right)

def add_to_gf_matrix(gf, tail, arr, nso, nmom, i_lev, j_lev, iwmax, num, n_occ, \
    zmsb, temp, left, right):
    """
    add_to_gf_matrix(gf, tail, arr, nso, nmom, i_lev, j_lev, iwmax, num, n_occ, \
        zmsb, temp, left, right)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 1188-1313
    
    Parameters
    ----------
    gf : complex array
    tail : complex array
    arr : int array
    nso : int
    nmom : int
    i_lev : int
    j_lev : int
    iwmax : int
    num : int
    n_occ : Occup_X100_Array
    	super-type
    
    zmsb : complex array
    temp : float
    left : bool
    right : bool
    
    """
    _hubbard_I.f90wrap_add_to_gf_matrix(gf=gf, tail=tail, arr=arr, nso=nso, \
        nmom=nmom, i_lev=i_lev, j_lev=j_lev, iwmax=iwmax, num=num, \
        n_occ=n_occ._handle, zmsb=zmsb, temp=temp, left=left, right=right)

def gf_hi_fullu(gf, tail, e0f, u, ummss, zmsb, nlm, iwmax, nmom, ns, temp, \
    verbosity, remove_split, n_lev):
    """
    atocc, atmag = gf_hi_fullu(gf, tail, e0f, u, ummss, zmsb, nlm, iwmax, nmom, ns, \
        temp, verbosity, remove_split, n_lev)
    
    
    Defined at \
        /mnt/beegfs/home/CPHT/leonid.poyurovskiy/DEVELOPMENT/MagInteract/development/MagInt_release/MagInt/fortran/MagInt/hubbard_I.f90 \
        lines 1315-1589
    
    Parameters
    ----------
    gf : complex array
    tail : complex array
    e0f : complex array
    u : float array
    ummss : float array
    zmsb : complex array
    nlm : int
    iwmax : int
    nmom : int
    ns : int
    temp : float
    verbosity : int
    remove_split : bool
    n_lev : int
    
    Returns
    -------
    atocc : float
    atmag : float
    
    """
    atocc, atmag = _hubbard_I.f90wrap_gf_hi_fullu(gf=gf, tail=tail, e0f=e0f, u=u, \
        ummss=ummss, zmsb=zmsb, nlm=nlm, iwmax=iwmax, nmom=nmom, ns=ns, temp=temp, \
        verbosity=verbosity, remove_split=remove_split, n_lev=n_lev)
    return atocc, atmag

