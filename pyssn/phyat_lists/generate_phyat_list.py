#!/usr/bin/env python

# Usage: python generate_phyat_list.py

import os
import sys
import subprocess
import argparse
from .manage_phyat_list import make_phyat_list, make_ionrec_file, merge_files, phyat2model
from .entries import execution_path
from ..fortran.runit import run_XSSN
import pyneb as pn

def config_pyneb():
    pn.atomicData.addDataFilePath(os.path.join(os.path.dirname(sys._getframe(1).f_code.co_filename), '../pyneb_data/'))
                                  
    # ***********************
    #
    # 3s1 mg2: -->s6 ok
    #
    # 3s2 mg1 si3: lk mg1, p4, cl_6(?) (out: 7 9 !) TBD**
    #
    #
    # 4s1 ca2  done incl. ca2 TBD**
    #pn.atomicData.setDataFile('k_i_atom_DP.dat')
    #pn.atomicData.setDataFile('k_i_coll_DP.dat')
    pn.atomicData.setDataFile('ca_ii_atom_DP.dat')
    pn.atomicData.setDataFile('ca_ii_coll_DP.dat')
    #
    # 4s2 ca1: TBD**
    #pn.atomicData.setDataFile('ca_i_atom_DP.dat')
    #pn.atomicData.setDataFile('ca_i_coll_DP.dat')
    #
    # 5s1 sr2: TBD** 
    # pn.atomicData.setDataFile('rb_i_atom_DP.dat')
    # pn.atomicData.setDataFile('rb_i_coll_DP.dat')
    #pn.atomicData.setDataFile('sr_ii_atom_DP.dat')
    #pn.atomicData.setDataFile('sr_ii_coll_DP.dat')
    #
    # 5s2 sr1: TBD**
    #pn.atomicData.setDataFile('sr_i_atom_DP.dat')
    #pn.atomicData.setDataFile('sr_i_coll_DP.dat')
    #
    # 5p6.6s1 ba2
    # (GS of La III = 5d elsewhere)
    # * ba_ii located after Fe_II ? chge GSCONFIG fe2 ?
    pn.atomicData.setDataFile('cs_i_atom_DP.dat')
    pn.atomicData.setDataFile('cs_i_coll_DP.dat')
    pn.atomicData.setDataFile('ba_ii_atom_DP.dat')
    pn.atomicData.setDataFile('ba_ii_coll_DP.dat')
    #
    # 6s2 ba1 (see again after doing ca1 and sr1)
    pn.atomicData.setDataFile('ba_i_atom_DP.dat')
    pn.atomicData.setDataFile('ba_i_coll_DP.dat')
    #
    # 2p1 ne6: lk 5 7 9 17 19 : done
    #
    pn.atomicData.setDataFile('f_v_atom_DP.dat')
    pn.atomicData.setDataFile('f_v_coll_DP.dat')
    pn.atomicData.setDataFile('na_vii_atom_DP.dat')
    pn.atomicData.setDataFile('na_vii_coll_DP.dat')
    pn.atomicData.setDataFile('al_ix_atom_DP.dat')
    pn.atomicData.setDataFile('al_ix_coll_DP.dat')
    pn.atomicData.setDataFile('sc_xvii_atom_DP.dat')
    pn.atomicData.setDataFile('sc_xvii_coll_DP.dat')
    pn.atomicData.setDataFile('v_xix_atom_DP.dat')
    pn.atomicData.setDataFile('v_xix_coll_DP.dat')
    #
    # 2p2 ne5 fe21: lk 4 16 18
    #
    #
    # 2p3 ne4 fe20: lk F3 Mg6? 15 17 19 20 !
    #
    #
    # 2p4 ne3 fe19: lk 2 7 14 16
    #
    #
    # 2p5 ne2, ar10, fe18: lk 1 34 6 10 11 13 15
    #
    #
    # 3p1 si2, ar6, fe14: lk 1 3 567 9 10 11 12 13 15 17 18(?) done incl. fe
    #
    # pn.atomicData.setDataFile('al_i_atom_DP.dat')
    # pn.atomicData.setDataFile('al_i_coll_DP.dat')
    # pn.atomicData.setDataFile('si_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('si_ii_coll_DP.dat')
    pn.atomicData.setDataFile('p_iii_atom_DP.dat')
    pn.atomicData.setDataFile('p_iii_coll_DP.dat')
    # pn.atomicData.setDataFile('s_iv_atom_DP.dat')
    # pn.atomicData.setDataFile('s_iv_coll_DP.dat')
    pn.atomicData.setDataFile('cl_v_atom_DP.dat')
    pn.atomicData.setDataFile('cl_v_coll_DP.dat')
    pn.atomicData.setDataFile('ar_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ar_vi_coll_DP.dat')
    pn.atomicData.setDataFile('k_vii_atom_DP.dat')
    pn.atomicData.setDataFile('k_vii_coll_DP.dat')
    pn.atomicData.setDataFile('ca_viii_atom_DP.dat')
    pn.atomicData.setDataFile('ca_viii_coll_DP.dat')
    pn.atomicData.setDataFile('sc_ix_atom_DP.dat')
    pn.atomicData.setDataFile('sc_ix_coll_DP.dat')
    pn.atomicData.setDataFile('ti_x_atom_DP.dat')
    pn.atomicData.setDataFile('ti_x_coll_DP.dat')
    pn.atomicData.setDataFile('v_xi_atom_DP.dat')
    pn.atomicData.setDataFile('v_xi_coll_DP.dat')
    pn.atomicData.setDataFile('cr_xii_atom_DP.dat')
    pn.atomicData.setDataFile('cr_xii_coll_DP.dat')
    pn.atomicData.setDataFile('mn_xiii_atom_DP.dat')
    pn.atomicData.setDataFile('mn_xiii_coll_DP.dat')
    pn.atomicData.setDataFile('fe_xiv_atom_DP.dat')
    pn.atomicData.setDataFile('fe_xiv_coll_DP.dat')
    pn.atomicData.setDataFile('co_xv_atom_DP.dat')
    pn.atomicData.setDataFile('co_xv_coll_DP.dat')
    pn.atomicData.setDataFile('ni_xvi_atom_DP.dat')
    pn.atomicData.setDataFile('ni_xvi_coll_DP.dat')
    pn.atomicData.setDataFile('cu_xvii_atom_DP.dat')
    pn.atomicData.setDataFile('cu_xvii_coll_DP.dat')
    pn.atomicData.setDataFile('zn_xviii_atom_DP.dat')
    pn.atomicData.setDataFile('zn_xviii_coll_DP.dat')
    #
    # 3p2 ar5, fe13: lk 1 2 789 10 11 12 14 16 17 TBD**
    #
    # *** 13 12 11 a remettre 
    # 3p5 ar2, fe10: lk 1 34567 11 13 14 ~ done incl. fe, xcpt zn_xiv bof & Omeg(cl_i)
    #
    pn.atomicData.setDataFile('cl_i_atom_DP.dat')
    pn.atomicData.setDataFile('cl_i_coll_DP.dat')
    pn.atomicData.setDataFile('ar_ii_atom_DP.dat')
    pn.atomicData.setDataFile('ar_ii_coll_DP.dat')
    pn.atomicData.setDataFile('k_iii_atom_DP.dat')
    pn.atomicData.setDataFile('k_iii_coll_DP.dat')
    pn.atomicData.setDataFile('ca_iv_atom_DP.dat')
    pn.atomicData.setDataFile('ca_iv_coll_DP.dat')
    pn.atomicData.setDataFile('sc_v_atom_DP.dat')
    pn.atomicData.setDataFile('sc_v_coll_DP.dat')
    pn.atomicData.setDataFile('ti_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ti_vi_coll_DP.dat')
    pn.atomicData.setDataFile('v_vii_atom_DP.dat')
    pn.atomicData.setDataFile('v_vii_coll_DP.dat')
    pn.atomicData.setDataFile('cr_viii_atom_DP.dat')
    pn.atomicData.setDataFile('cr_viii_coll_DP.dat')
    pn.atomicData.setDataFile('mn_ix_atom_DP.dat')
    pn.atomicData.setDataFile('mn_ix_coll_DP.dat')
    pn.atomicData.setDataFile('fe_x_atom_DP.dat')
    pn.atomicData.setDataFile('fe_x_coll_DP.dat')
    pn.atomicData.setDataFile('co_xi_atom_DP.dat')
    pn.atomicData.setDataFile('co_xi_coll_DP.dat')
    pn.atomicData.setDataFile('ni_xii_atom_DP.dat')
    pn.atomicData.setDataFile('ni_xii_coll_DP.dat')
    pn.atomicData.setDataFile('cu_xiii_atom_DP.dat')
    pn.atomicData.setDataFile('cu_xiii_coll_DP.dat')
    pn.atomicData.setDataFile('zn_xiv_atom_DP.dat')
    pn.atomicData.setDataFile('zn_xiv_coll_DP.dat')
    # 
    # 4p1 kr6:
    #
    pn.atomicData.setDataFile('ga_i_atom_DP.dat')
    pn.atomicData.setDataFile('ga_i_coll_DP.dat')
    pn.atomicData.setDataFile('ge_ii_atom_DP.dat')
    pn.atomicData.setDataFile('ge_ii_coll_DP.dat')
    pn.atomicData.setDataFile('as_iii_atom_DP.dat')
    pn.atomicData.setDataFile('as_iii_coll_DP.dat')
    pn.atomicData.setDataFile('se_iv_atom_DP.dat')
    pn.atomicData.setDataFile('se_iv_coll_DP.dat')
    pn.atomicData.setDataFile('br_v_atom_DP.dat')
    pn.atomicData.setDataFile('br_v_coll_DP.dat')
    pn.atomicData.setDataFile('kr_vi_atom_DP.dat')
    pn.atomicData.setDataFile('kr_vi_coll_DP.dat')
    pn.atomicData.setDataFile('rb_vii_atom_DP.dat')
    pn.atomicData.setDataFile('rb_vii_coll_DP.dat')
    pn.atomicData.setDataFile('sr_viii_atom_DP.dat')
    pn.atomicData.setDataFile('sr_viii_coll_DP.dat')
    pn.atomicData.setDataFile('y_ix_atom_DP.dat')
    pn.atomicData.setDataFile('y_ix_coll_DP.dat')
    pn.atomicData.setDataFile('zr_x_atom_DP.dat')
    pn.atomicData.setDataFile('zr_x_coll_DP.dat')
    pn.atomicData.setDataFile('nb_xi_atom_DP.dat')
    pn.atomicData.setDataFile('nb_xi_coll_DP.dat')
    pn.atomicData.setDataFile('mo_xii_atom_DP.dat')
    pn.atomicData.setDataFile('mo_xii_coll_DP.dat')
    #
    # 4p2 kr5: lk 1 2 3pb 4pb 5only_redone 6pbcoll 7...
    #
    pn.atomicData.setDataFile('ge_i_atom_DP.dat')
    pn.atomicData.setDataFile('ge_i_coll_DP.dat')
    pn.atomicData.setDataFile('as_ii_atom_DP.dat')
    pn.atomicData.setDataFile('as_ii_coll_DP.dat')
    pn.atomicData.setDataFile('se_iii_atom_DP.dat')
    pn.atomicData.setDataFile('se_iii_coll_DP.dat')
    pn.atomicData.setDataFile('br_iv_atom_DP.dat')
    pn.atomicData.setDataFile('br_iv_coll_DP.dat')
    pn.atomicData.setDataFile('kr_v_atom_DP.dat')
    pn.atomicData.setDataFile('kr_v_coll_DP.dat')
    pn.atomicData.setDataFile('rb_vi_atom_DP.dat')
    pn.atomicData.setDataFile('rb_vi_coll_DP.dat')
    pn.atomicData.setDataFile('sr_vii_atom_DP.dat')
    pn.atomicData.setDataFile('sr_vii_coll_DP.dat')
    pn.atomicData.setDataFile('y_viii_atom_DP.dat')
    pn.atomicData.setDataFile('y_viii_coll_DP.dat')
    pn.atomicData.setDataFile('zr_ix_atom_DP.dat')
    pn.atomicData.setDataFile('zr_ix_coll_DP.dat')
    pn.atomicData.setDataFile('nb_x_atom_DP.dat')
    pn.atomicData.setDataFile('nb_x_coll_DP.dat')
    pn.atomicData.setDataFile('mo_xi_atom_DP.dat')
    pn.atomicData.setDataFile('mo_xi_coll_DP.dat')
    #
    # 4p3 kr4: lk 1 2 345_pb_scalOm_decimAij 6..
    #
    pn.atomicData.setDataFile('as_i_atom_DP.dat')
    pn.atomicData.setDataFile('as_i_coll_DP.dat')
    pn.atomicData.setDataFile('se_ii_atom_DP.dat')
    pn.atomicData.setDataFile('se_ii_coll_DP.dat')
    pn.atomicData.setDataFile('br_iii_atom_DP.dat')
    pn.atomicData.setDataFile('br_iii_coll_DP.dat')
    pn.atomicData.setDataFile('kr_iv_atom_DP.dat')
    pn.atomicData.setDataFile('kr_iv_coll_DP.dat')
    pn.atomicData.setDataFile('rb_v_atom_DP.dat')
    pn.atomicData.setDataFile('rb_v_coll_DP.dat')
    pn.atomicData.setDataFile('sr_vi_atom_DP.dat')
    pn.atomicData.setDataFile('sr_vi_coll_DP.dat')
    pn.atomicData.setDataFile('y_vii_atom_DP.dat')
    pn.atomicData.setDataFile('y_vii_coll_DP.dat')
    pn.atomicData.setDataFile('zr_viii_atom_DP.dat')
    pn.atomicData.setDataFile('zr_viii_coll_DP.dat')
    pn.atomicData.setDataFile('nb_ix_atom_DP.dat')
    pn.atomicData.setDataFile('nb_ix_coll_DP.dat')
    pn.atomicData.setDataFile('mo_x_atom_DP.dat')
    pn.atomicData.setDataFile('mo_x_coll_DP.dat')
    #
    # 4p4 kr3:  lk 1 2 kr3_decimAij rb4pb_scal 5 6..
    #
    pn.atomicData.setDataFile('se_i_atom_DP.dat')
    pn.atomicData.setDataFile('se_i_coll_DP.dat')
    pn.atomicData.setDataFile('br_ii_atom_DP.dat')
    pn.atomicData.setDataFile('br_ii_coll_DP.dat')
    pn.atomicData.setDataFile('kr_iii_atom_DP.dat')
    pn.atomicData.setDataFile('kr_iii_coll_DP.dat')
    pn.atomicData.setDataFile('rb_iv_atom_DP.dat')
    pn.atomicData.setDataFile('rb_iv_coll_DP.dat')
    pn.atomicData.setDataFile('sr_v_atom_DP.dat')
    pn.atomicData.setDataFile('sr_v_coll_DP.dat')
    pn.atomicData.setDataFile('y_vi_atom_DP.dat')
    pn.atomicData.setDataFile('y_vi_coll_DP.dat')
    pn.atomicData.setDataFile('zr_vii_atom_DP.dat')
    pn.atomicData.setDataFile('zr_vii_coll_DP.dat')
    pn.atomicData.setDataFile('nb_viii_atom_DP.dat')
    pn.atomicData.setDataFile('nb_viii_coll_DP.dat')
    pn.atomicData.setDataFile('mo_ix_atom_DP.dat')
    pn.atomicData.setDataFile('mo_ix_coll_DP.dat')
    #
    # 4p5 kr2:
    #
    pn.atomicData.setDataFile('br_i_atom_DP.dat')
    pn.atomicData.setDataFile('br_i_coll_DP.dat')
    pn.atomicData.setDataFile('kr_ii_atom_DP.dat')
    pn.atomicData.setDataFile('kr_ii_coll_DP.dat')
    pn.atomicData.setDataFile('rb_iii_atom_DP.dat')
    pn.atomicData.setDataFile('rb_iii_coll_DP.dat')
    pn.atomicData.setDataFile('sr_iv_atom_DP.dat')
    pn.atomicData.setDataFile('sr_iv_coll_DP.dat')
    pn.atomicData.setDataFile('y_v_atom_DP.dat')
    pn.atomicData.setDataFile('y_v_coll_DP.dat')
    pn.atomicData.setDataFile('zr_vi_atom_DP.dat')
    pn.atomicData.setDataFile('zr_vi_coll_DP.dat')
    pn.atomicData.setDataFile('nb_vii_atom_DP.dat')
    pn.atomicData.setDataFile('nb_vii_coll_DP.dat')
    pn.atomicData.setDataFile('mo_viii_atom_DP.dat')
    pn.atomicData.setDataFile('mo_viii_coll_DP.dat')
    #
    # 5p1 xe6: lk all but xe6
    #
    pn.atomicData.setDataFile('in_i_atom_DP.dat')
    pn.atomicData.setDataFile('in_i_coll_DP.dat')
    pn.atomicData.setDataFile('sn_ii_atom_DP.dat')
    pn.atomicData.setDataFile('sn_ii_coll_DP.dat')
    pn.atomicData.setDataFile('sb_iii_atom_DP.dat')
    pn.atomicData.setDataFile('sb_iii_coll_DP.dat')
    pn.atomicData.setDataFile('te_iv_atom_DP.dat')
    pn.atomicData.setDataFile('te_iv_coll_DP.dat')
    pn.atomicData.setDataFile('i_v_atom_DP.dat')
    pn.atomicData.setDataFile('i_v_coll_DP.dat')
    pn.atomicData.setDataFile('xe_vi_atom_DP.dat')
    pn.atomicData.setDataFile('xe_vi_coll_DP.dat')
    pn.atomicData.setDataFile('cs_vii_atom_DP.dat')
    pn.atomicData.setDataFile('cs_vii_coll_DP.dat')
    pn.atomicData.setDataFile('ba_viii_atom_DP.dat')
    pn.atomicData.setDataFile('ba_viii_coll_DP.dat')
    pn.atomicData.setDataFile('la_ix_atom_DP.dat')
    pn.atomicData.setDataFile('la_ix_coll_DP.dat')
    #
    # 5p2 xe5:
    #
    pn.atomicData.setDataFile('sn_i_atom_DP.dat')
    pn.atomicData.setDataFile('sn_i_coll_DP.dat')
    pn.atomicData.setDataFile('sb_ii_atom_DP.dat')
    pn.atomicData.setDataFile('sb_ii_coll_DP.dat')
    pn.atomicData.setDataFile('te_iii_atom_DP.dat')
    pn.atomicData.setDataFile('te_iii_coll_DP.dat')
    pn.atomicData.setDataFile('i_iv_atom_DP.dat')
    pn.atomicData.setDataFile('i_iv_coll_DP.dat')
    pn.atomicData.setDataFile('xe_v_atom_DP.dat')
    pn.atomicData.setDataFile('xe_v_coll_DP.dat')
    pn.atomicData.setDataFile('cs_vi_atom_DP.dat')
    pn.atomicData.setDataFile('cs_vi_coll_DP.dat')
    pn.atomicData.setDataFile('ba_vii_atom_DP.dat')
    pn.atomicData.setDataFile('ba_vii_coll_DP.dat')
    pn.atomicData.setDataFile('la_viii_atom_DP.dat')
    pn.atomicData.setDataFile('la_viii_coll_DP.dat')
    pn.atomicData.setDataFile('ce_ix_atom_DP.dat')
    pn.atomicData.setDataFile('ce_ix_coll_DP.dat')
    #
    # 5p3 xe4:  lk all but xe4_decAij_reconsidered
    #
    pn.atomicData.setDataFile('sb_i_atom_DP.dat')
    pn.atomicData.setDataFile('sb_i_coll_DP.dat')
    pn.atomicData.setDataFile('te_ii_atom_DP.dat')
    pn.atomicData.setDataFile('te_ii_coll_DP.dat')
    pn.atomicData.setDataFile('i_iii_atom_DP.dat')
    pn.atomicData.setDataFile('i_iii_coll_DP.dat')
    pn.atomicData.setDataFile('xe_iv_atom_DP.dat')
    pn.atomicData.setDataFile('xe_iv_coll_DP.dat')
    pn.atomicData.setDataFile('cs_v_atom_DP.dat')
    pn.atomicData.setDataFile('cs_v_coll_DP.dat')
    pn.atomicData.setDataFile('ba_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ba_vi_coll_DP.dat')
    pn.atomicData.setDataFile('la_vii_atom_DP.dat')
    pn.atomicData.setDataFile('la_vii_coll_DP.dat')
    pn.atomicData.setDataFile('ce_viii_atom_DP.dat')
    pn.atomicData.setDataFile('ce_viii_coll_DP.dat')
    #
    # 5p4 xe3:  lk all but xe3_reconsidered 
    #
    pn.atomicData.setDataFile('te_i_atom_DP.dat')
    pn.atomicData.setDataFile('te_i_coll_DP.dat')
    pn.atomicData.setDataFile('i_ii_atom_DP.dat')
    pn.atomicData.setDataFile('i_ii_coll_DP.dat')
    pn.atomicData.setDataFile('xe_iii_atom_DP.dat')
    pn.atomicData.setDataFile('xe_iii_coll_DP.dat')
    pn.atomicData.setDataFile('cs_iv_atom_DP.dat')
    pn.atomicData.setDataFile('cs_iv_coll_DP.dat')
    pn.atomicData.setDataFile('ba_v_atom_DP.dat')
    pn.atomicData.setDataFile('ba_v_coll_DP.dat')
    pn.atomicData.setDataFile('la_vi_atom_DP.dat')
    pn.atomicData.setDataFile('la_vi_coll_DP.dat')
    pn.atomicData.setDataFile('ce_vii_atom_DP.dat')
    pn.atomicData.setDataFile('ce_vii_coll_DP.dat')
    #
    # 5p5 xe2:
    #
    pn.atomicData.setDataFile('i_i_atom_DP.dat')
    pn.atomicData.setDataFile('i_i_coll_DP.dat')
    pn.atomicData.setDataFile('xe_ii_atom_DP.dat')
    pn.atomicData.setDataFile('xe_ii_coll_DP.dat')
    pn.atomicData.setDataFile('cs_iii_atom_DP.dat')
    pn.atomicData.setDataFile('cs_iii_coll_DP.dat')
    pn.atomicData.setDataFile('ba_iv_atom_DP.dat')
    pn.atomicData.setDataFile('ba_iv_coll_DP.dat')
    pn.atomicData.setDataFile('la_v_atom_DP.dat')
    pn.atomicData.setDataFile('la_v_coll_DP.dat')
    pn.atomicData.setDataFile('ce_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ce_vi_coll_DP.dat')
    #
    # 6p1 pb2:
    #
    pn.atomicData.setDataFile('tl_i_atom_DP.dat')
    pn.atomicData.setDataFile('tl_i_coll_DP.dat')
    pn.atomicData.setDataFile('pb_ii_atom_DP.dat')
    pn.atomicData.setDataFile('pb_ii_coll_DP.dat')
    pn.atomicData.setDataFile('bi_iii_atom_DP.dat')
    pn.atomicData.setDataFile('bi_iii_coll_DP.dat')
    #
    # 6p2 bi2:
    #
    pn.atomicData.setDataFile('pb_i_atom_DP.dat')
    pn.atomicData.setDataFile('pb_i_coll_DP.dat')
    pn.atomicData.setDataFile('bi_ii_atom_DP.dat')
    pn.atomicData.setDataFile('bi_ii_coll_DP.dat')
    # 
    # 3d1: fe_viii lk all but fe. done incl. fe
    #
    pn.atomicData.setDataFile('sc_iii_atom_DP.dat')
    pn.atomicData.setDataFile('sc_iii_coll_DP.dat')
    pn.atomicData.setDataFile('ti_iv_atom_DP.dat')
    pn.atomicData.setDataFile('ti_iv_coll_DP.dat')
    pn.atomicData.setDataFile('v_v_atom_DP.dat')
    pn.atomicData.setDataFile('v_v_coll_DP.dat')
    pn.atomicData.setDataFile('cr_vi_atom_DP.dat')
    pn.atomicData.setDataFile('cr_vi_coll_DP.dat')
    pn.atomicData.setDataFile('mn_vii_atom_DP.dat')
    pn.atomicData.setDataFile('mn_vii_coll_DP.dat')
    pn.atomicData.setDataFile('fe_viii_atom_DP.dat')
    pn.atomicData.setDataFile('fe_viii_coll_DP.dat')
    pn.atomicData.setDataFile('co_ix_atom_DP.dat')
    pn.atomicData.setDataFile('co_ix_coll_DP.dat')
    pn.atomicData.setDataFile('ni_x_atom_DP.dat')
    pn.atomicData.setDataFile('ni_x_coll_DP.dat')
    pn.atomicData.setDataFile('cu_xi_atom_DP.dat')
    pn.atomicData.setDataFile('cu_xi_coll_DP.dat')
    pn.atomicData.setDataFile('zn_xii_atom_DP.dat')
    pn.atomicData.setDataFile('zn_xii_coll_DP.dat')
    #
    # 3d2: fe_vii lk all but fe. done
    #  ***28/7/17 sc_ii abst ? co_viii corr ressayer ?
    pn.atomicData.setDataFile('sc_ii_atom_DP.dat')
    pn.atomicData.setDataFile('sc_ii_coll_DP.dat')
    pn.atomicData.setDataFile('ti_iii_atom_DP.dat')
    pn.atomicData.setDataFile('ti_iii_coll_DP.dat')
    pn.atomicData.setDataFile('v_iv_atom_DP.dat')
    pn.atomicData.setDataFile('v_iv_coll_DP.dat')
    pn.atomicData.setDataFile('cr_v_atom_DP.dat')
    pn.atomicData.setDataFile('cr_v_coll_DP.dat')
    pn.atomicData.setDataFile('mn_vi_atom_DP.dat')
    pn.atomicData.setDataFile('mn_vi_coll_DP.dat')
    # pn.atomicData.setDataFile('fe_vii_atom_WB08.dat')
    # pn.atomicData.setDataFile('fe_vii_coll_WB08.dat')
    pn.atomicData.setDataFile('co_viii_atom_DP.dat')
    pn.atomicData.setDataFile('co_viii_coll_DP.dat')
    pn.atomicData.setDataFile('ni_ix_atom_DP.dat')
    pn.atomicData.setDataFile('ni_ix_coll_DP.dat')
    pn.atomicData.setDataFile('cu_x_atom_DP.dat')
    pn.atomicData.setDataFile('cu_x_coll_DP.dat')
    pn.atomicData.setDataFile('zn_xi_atom_DP.dat')
    pn.atomicData.setDataFile('zn_xi_coll_DP.dat')
    #
    # 3d3: fe_vi lk all but fe. done xcpt ti_ii
    #
    # pn.atomicData.setDataFile('ti_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('ti_ii_coll_DP.dat')
    pn.atomicData.setDataFile('v_iii_atom_DP.dat')
    pn.atomicData.setDataFile('v_iii_coll_DP.dat')
    pn.atomicData.setDataFile('cr_iv_atom_DP.dat')
    pn.atomicData.setDataFile('cr_iv_coll_DP.dat')
    pn.atomicData.setDataFile('mn_v_atom_DP.dat')
    pn.atomicData.setDataFile('mn_v_coll_DP.dat')
    # pn.atomicData.setDataFile('fe_vi')
    # pn.atomicData.setDataFile('fe_vi')
    pn.atomicData.setDataFile('co_vii_atom_DP.dat')
    pn.atomicData.setDataFile('co_vii_coll_DP.dat')
    pn.atomicData.setDataFile('ni_viii_atom_DP.dat')
    pn.atomicData.setDataFile('ni_viii_coll_DP.dat')
    #
    # 3d4: fe_v lk all but fe. done incl. fe, xcpt v_ii
    #
    # pn.atomicData.setDataFile('v_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('v_ii_coll_DP.dat')
    pn.atomicData.setDataFile('cr_iii_atom_DP.dat')
    pn.atomicData.setDataFile('cr_iii_coll_DP.dat')
    pn.atomicData.setDataFile('mn_iv_atom_DP.dat')
    pn.atomicData.setDataFile('mn_iv_coll_DP.dat')
    pn.atomicData.setDataFile('fe_v_atom_DP.dat')
    pn.atomicData.setDataFile('fe_v_coll_DP.dat')
    pn.atomicData.setDataFile('co_vi_atom_DP.dat')
    pn.atomicData.setDataFile('co_vi_coll_DP.dat')
    pn.atomicData.setDataFile('ni_vii_atom_DP.dat')
    pn.atomicData.setDataFile('ni_vii_coll_DP.dat')
    #
    # 3d5: fe_iv lk all but fe. done incl. fe
    #  **cr_ii in the making ?
    # pn.atomicData.setDataFile('cr_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('cr_ii_coll_DP.dat')
    pn.atomicData.setDataFile('mn_iii_atom_DP.dat')
    pn.atomicData.setDataFile('mn_iii_coll_DP.dat')
    pn.atomicData.setDataFile('fe_iv_atom_DP.dat')
    pn.atomicData.setDataFile('fe_iv_coll_DP.dat')
    pn.atomicData.setDataFile('co_v_atom_DP.dat')
    pn.atomicData.setDataFile('co_v_coll_DP.dat')
    pn.atomicData.setDataFile('ni_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ni_vi_coll_DP.dat')
    #
    # 3d6: fe_iii lk all but fe. done incl. fe
    #  **? mn_ii (4s ?) entre ca_ii et fe_ii
    #  **? zn_vii abst
    pn.atomicData.setDataFile('mn_ii_atom_DP.dat')
    pn.atomicData.setDataFile('mn_ii_coll_DP.dat')
    # pn.atomicData.setDataFile('fe_iii_atom_....dat')
    # pn.atomicData.setDataFile('fe_iii_coll_....dat')
    pn.atomicData.setDataFile('co_iv_atom_DP.dat')
    pn.atomicData.setDataFile('co_iv_coll_DP.dat')
    pn.atomicData.setDataFile('ni_v_atom_DP.dat')
    pn.atomicData.setDataFile('ni_v_coll_DP.dat')
    pn.atomicData.setDataFile('cu_vi_atom_DP.dat')
    pn.atomicData.setDataFile('cu_vi_coll_DP.dat')
    pn.atomicData.setDataFile('zn_vii_atom_DP.dat')
    pn.atomicData.setDataFile('zn_vii_coll_DP.dat')
    #
    # 3d7: (fe_ii) ni_iv  lk all but fe. done xcpt ga_vii ge_viii
    #  ** fe_ii elsewhere (4s)
    #  ** ni_iv cu_v abst ? 
    # pn.atomicData.setDataFile('fe_ii_atom_....dat')
    # pn.atomicData.setDataFile('fe_ii_coll_....dat')
    pn.atomicData.setDataFile('co_iii_atom_DP.dat')
    pn.atomicData.setDataFile('co_iii_coll_DP.dat')
    pn.atomicData.setDataFile('ni_iv_atom_DP.dat')
    pn.atomicData.setDataFile('ni_iv_coll_DP.dat')
    #pn.atomicData.setDataFile('cu_v_atom_DP.dat')
    #pn.atomicData.setDataFile('cu_v_coll_DP.dat')
    pn.atomicData.setDataFile('zn_vi_atom_DP.dat')
    pn.atomicData.setDataFile('zn_vi_coll_DP.dat')
    # pn.atomicData.setDataFile('ga_vii_atom_DP.dat')
    pn.atomicData.setDataFile('ga_vii_coll_DP.dat')
    # pn.atomicData.setDataFile('ge_viii_atom_DP.dat')
    pn.atomicData.setDataFile('ge_viii_coll_DP.dat')
    #
    # 3d8: (fe_i) ni_iii  lk all but ni. done incl. ni, xcpt ga_vi ge_vii
    # ** co_ii abst ?
    # ** a finir : ga_vi ge_vii
    pn.atomicData.setDataFile('co_ii_atom_DP.dat')
    pn.atomicData.setDataFile('co_ii_coll_DP.dat')
    pn.atomicData.setDataFile('ni_iii_atom_DP.dat')
    pn.atomicData.setDataFile('ni_iii_coll_DP.dat')
    pn.atomicData.setDataFile('cu_iv_atom_DP.dat')
    pn.atomicData.setDataFile('cu_iv_coll_DP.dat')
    pn.atomicData.setDataFile('zn_v_atom_DP.dat')
    pn.atomicData.setDataFile('zn_v_coll_DP.dat')
    # pn.atomicData.setDataFile('ga_vi_atom_DP.dat')
    # pn.atomicData.setDataFile('ga_vi_coll_DP.dat')
    # pn.atomicData.setDataFile('ge_vii_atom_DP.dat')
    # pn.atomicData.setDataFile('ge_vii_coll_DP.dat')
    #
    # 3d9: ni_ii, zn_iv  lk all but ni. done incl. ni
    # ** ni_ii abst (seul non 2 niv ?),  sr_xii corrige a verif
    pn.atomicData.setDataFile('ni_ii_atom_DP.dat') 
    pn.atomicData.setDataFile('ni_ii_coll_DP.dat')
    pn.atomicData.setDataFile('cu_iii_atom_DP.dat')
    pn.atomicData.setDataFile('cu_iii_coll_DP.dat')
    pn.atomicData.setDataFile('zn_iv_atom_DP.dat')
    pn.atomicData.setDataFile('zn_iv_coll_DP.dat')
    pn.atomicData.setDataFile('ga_v_atom_DP.dat')
    pn.atomicData.setDataFile('ga_v_coll_DP.dat')
    pn.atomicData.setDataFile('ge_vi_atom_DP.dat')
    pn.atomicData.setDataFile('ge_vi_coll_DP.dat')
    pn.atomicData.setDataFile('as_vii_atom_DP.dat')
    pn.atomicData.setDataFile('as_vii_coll_DP.dat')
    pn.atomicData.setDataFile('se_viii_atom_DP.dat')
    pn.atomicData.setDataFile('se_viii_coll_DP.dat')
    pn.atomicData.setDataFile('br_ix_atom_DP.dat')
    pn.atomicData.setDataFile('br_ix_coll_DP.dat')
    pn.atomicData.setDataFile('kr_x_atom_DP.dat')
    pn.atomicData.setDataFile('kr_x_coll_DP.dat')
    pn.atomicData.setDataFile('rb_xi_atom_DP.dat')
    pn.atomicData.setDataFile('rb_xi_coll_DP.dat')
    pn.atomicData.setDataFile('sr_xii_atom_DP.dat')
    pn.atomicData.setDataFile('sr_xii_coll_DP.dat')
    #
    # 3d10: cu_ii
    # (**Rem: In_V juste avt CuII ?)
    pn.atomicData.setDataFile('cu_ii_atom_DP.dat')
    pn.atomicData.setDataFile('cu_ii_coll_DP.dat')
    pn.atomicData.setDataFile('zn_iii_atom_DP.dat')
    pn.atomicData.setDataFile('zn_iii_coll_DP.dat')
    #
    # 4d1: y_iii
    #
    #
    # 4d7: xe_xii
    #
    #
    # 4d8: xe_xi
    #
    # pn.atomicData.setDataFile('rh_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('rh_ii_coll_DP.dat')
    # pn.atomicData.setDataFile('pd_iii_atom_DP.dat')
    # pn.atomicData.setDataFile('pd_iii_coll_DP.dat')
    # pn.atomicData.setDataFile('ag_iv_atom_DP.dat')
    # pn.atomicData.setDataFile('ag_iv_coll_DP.dat')
    """
    pn.atomicData.setDataFile('cd_v_atom_DP.dat')
    pn.atomicData.setDataFile('cd_v_coll_DP.dat')
    pn.atomicData.setDataFile('in_vi_atom_DP.dat')
    pn.atomicData.setDataFile('in_vi_coll_DP.dat')
    pn.atomicData.setDataFile('sn_vii_atom_DP.dat')
    pn.atomicData.setDataFile('sn_vii_coll_DP.dat')
    pn.atomicData.setDataFile('sb_viii_atom_DP.dat')
    pn.atomicData.setDataFile('sb_viii_coll_DP.dat')
    pn.atomicData.setDataFile('te_ix_atom_DP.dat')
    pn.atomicData.setDataFile('te_ix_coll_DP.dat')
    pn.atomicData.setDataFile('i_x_atom_DP.dat')
    pn.atomicData.setDataFile('i_x_coll_DP.dat')
    pn.atomicData.setDataFile('xe_xi_atom_DP.dat')
    pn.atomicData.setDataFile('xe_xi_coll_DP.dat')
    pn.atomicData.setDataFile('cs_xii_atom_DP.dat')
    pn.atomicData.setDataFile('cs_xii_coll_DP.dat')
    pn.atomicData.setDataFile('ba_xiii_atom_DP.dat')
    pn.atomicData.setDataFile('ba_xiii_coll_DP.dat')
    """
    #
    # 4d9: xe_x
    #
    # pn.atomicData.setDataFile('pd_ii_atom_DP.dat')
    # pn.atomicData.setDataFile('pd_ii_coll_DP.dat')
    # pn.atomicData.setDataFile('ag_iii_atom_DP.dat')
    # pn.atomicData.setDataFile('ag_iii_coll_DP.dat')
    pn.atomicData.setDataFile('cd_iv_atom_DP.dat')
    pn.atomicData.setDataFile('cd_iv_coll_DP.dat')
    pn.atomicData.setDataFile('in_v_atom_DP.dat')
    pn.atomicData.setDataFile('in_v_coll_DP.dat')
    #pn.atomicData.setDataFile('sn_vi_atom_DP.dat')
    #pn.atomicData.setDataFile('sn_vi_coll_DP.dat')
    pn.atomicData.setDataFile('sb_vii_atom_DP.dat')
    pn.atomicData.setDataFile('sb_vii_coll_DP.dat')
    pn.atomicData.setDataFile('te_viii_atom_DP.dat')
    pn.atomicData.setDataFile('te_viii_coll_DP.dat')
    pn.atomicData.setDataFile('i_ix_atom_DP.dat')
    pn.atomicData.setDataFile('i_ix_coll_DP.dat')
    pn.atomicData.setDataFile('xe_x_atom_DP.dat')
    pn.atomicData.setDataFile('xe_x_coll_DP.dat')
    pn.atomicData.setDataFile('cs_xi_atom_DP.dat')
    pn.atomicData.setDataFile('cs_xi_coll_DP.dat')
    pn.atomicData.setDataFile('ba_xii_atom_DP.dat')
    pn.atomicData.setDataFile('ba_xii_coll_DP.dat')
    #
    # 5p6.4f ce_iv (sf la_iii: 5d1)
    #
    pn.atomicData.setDataFile('la_iii_atom_DP.dat')
    pn.atomicData.setDataFile('la_iii_coll_DP.dat')
    pn.atomicData.setDataFile('ce_iv_atom_DP.dat')
    pn.atomicData.setDataFile('ce_iv_coll_DP.dat')
    pn.atomicData.setDataFile('pr_v_atom_DP.dat')
    pn.atomicData.setDataFile('pr_v_coll_DP.dat')
    pn.atomicData.setDataFile('nd_vi_atom_DP.dat')
    pn.atomicData.setDataFile('nd_vi_coll_DP.dat')
    # 



# ***********************
# Set to None (atoms = None) to generate all of the available atoms
atoms = None
# atoms = ['C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2', 'Fe3']

# atoms = ['F5', 'Ne6', 'Na7', 'Al9', 'Sc17', 'V19', 'Ga1', 'Ge1', 'Ge2', 'As1', 'As2', 'As3', 'Se1', 'Se2', 'Se3', 'Se4', 'Br1', 'Br2', 'Br3', 'Br4', 'Br5', 'Kr2', 'Kr3', 'Kr4', 'Kr5', 'Kr6', 'Rb3', 'Rb4', 'Rb5', 'Rb6', 'Rb7', 'Sr4', 'Sr5', 'Sr6', 'Sr7', 'Sr8', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9', 'Zr6', 'Zr7', 'Zr8', 'Zr9', 'Zr10', 'Nb7', 'Nb8', 'Nb9', 'Nb10', 'Nb11', 'Mo8', 'Mo9', 'Mo10', 'Mo11', 'Mo12','In1','Sn1', 'Sn2', 'Sb1', 'Sb2', 'Sb3', 'Te1', 'Te2', 'Te3', 'Te4', 'I1', 'I2', 'I3', 'I4', 'I5', 'Xe2', 'Xe3', 'Xe4', 'Xe5', 'Xe6', 'Cs3', 'Cs4', 'Cs5', 'Cs6', 'Cs7', 'Ba4', 'Ba5', 'Ba6', 'Ba7', 'Ba8', 'La5', 'L6', 'La7', 'La8', 'La9', 'Ce6', 'Ce7', 'Ce8', 'Ce9', 'Tl1', 'Pb2', 'Bi3', 'Pb1', 'Bi2', 'Ba2']

# File containing additional data in the phyat_list format
extra_file = None
# extra_file = 'phyat_list_DP_p3.dat'
# extra_file = 'phyat_list_DP_nouv.dat'
# extra_file = 'phyat_list_bidon.dat'
# extra_file = 'phyat_list_DP_nouv.dat'

# Name of the output
# filename = 'listep3t1_12.dat'
# filename = 'listep3SIIt4_10.dat'
# filename = 'listep3SIIt1_8.dat'
# filename = 'listeCIIIt1.5_4.6.4.dat'
# filename = 'liste_test_DP.dat'
# filename = 'liste_extrafile_none_noBaIIgenerate.dat'
filename = 'liste_phyat_coll.dat'
phy_cond_file = 'phy_cond.dat'

ref_lines_dic = {}

NLevels_dic = {}

up_lev_rule_dic = {}

Del_ion = []

# Transitions set to 0.0
Aij_zero_dic = {}


def make_all_lists():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--abund_file", help="Abundances file")
    parser.add_argument("-I", "--ion_frac_file", help="Ionic fractions file")
    parser.add_argument("-P", "--phy_cond_file", help="Physical conditions file")
    parser.add_argument("-C", "--outputcond_file", help="Output conditions file")
    parser.add_argument("-O", "--phyat_file", help="Output phyat file")
    
    args = parser.parse_args()
    
    make_ionrec_file(abund_file=args.abund_file, ion_frac_file=args.ion_frac_file, out_file=execution_path('ions_rec.dat', extra='../fortran/'))
    print('ionrec file done')
    run_XSSN(outputcond_file = args.outputcond_file)
    print('XSSN run')
    config_pyneb()
    print('PyNeb configured')
    make_phyat_list(filename=execution_path('liste_phyat_coll.dat'), atoms=None, cut=1e-4, E_cut=20,
          verbose=False, notry=False, NLevels=50, 
          ref_lines_dic=ref_lines_dic,
          NLevels_dic=NLevels_dic,
          up_lev_rule_dic=up_lev_rule_dic,
          Aij_zero_dic=Aij_zero_dic,
          Del_ion = Del_ion,
          phy_cond_file = args.phy_cond_file,
          extra_file=extra_file)
    print('phyat col done')
    
    merge_files((execution_path('liste_phyat_rec.dat', extra='../fortran/'), execution_path('liste_phyat_coll.dat'), 
                 execution_path('liste_phyat_others.dat')), args.phyat_file)
    print('Files merged')
    
def make_list_model():
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--abund_file", help="Abundances file")
    parser.add_argument("-I", "--ion_frac_file", help="Ionic fractions file")
    parser.add_argument("-O", "--phyat_file", help="Output phyat file")
    parser.add_argument("-M", "--model_file", help="Output model file")
    parser.add_argument("-N", "--norm_hbeta", help="Hbeta value", default=10000)
    parser.add_argument("-F", "--ion_frac_min", help="Ion Frac min", default=0.0001)
    
    args = parser.parse_args()
    
        
    phyat2model(args.phyat_file, args.model_file, norm_hbeta=float(args.norm_hbeta), ion_frac_file=args.ion_frac_file, 
                abund_file=args.abund_file, ion_frac_min=float(args.ion_frac_min))
    print('Model done')
    
    
    
"""
 --abund_file=asplund_2009.dat --ion_frac_file=50_-2_R_ionfrac.dat --phy_cond_file=phy_cond.dat --outputcond_file=outputcond.dat` 
`--phyat_file=liste_phyat_test1.dat --model_file=liste_model_test1.dat --norm_hbeta=10000 --ion_frac_min=0.0001`
"""
