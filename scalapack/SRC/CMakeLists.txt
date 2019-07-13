set (ALLAUX  pjlaenv.f  pilaenvx.f piparmq.f pilaver.f pmpim2.f pmpcol.f
     ilaenv.f ilaenv2stage.f ieeeck.f lsamen.f iparmq.f iparam2stage.F
     ilaprec.f ilatrans.f ilauplo.f iladiag.f chla_transtype.f
     ../INSTALL/ilaver.f ../INSTALL/lsame.f xerbla.f xerbla_array.f
     ../INSTALL/slamch.f)

set (ALLAUX-C  pbchkvect.c getpbbuf.c pcrot.c pslaiect.c pdlaiect.c pzrot.c slamov.c clamov.c dlamov.c zlamov.c)

set (DZLAUX  
   pdlabad.f pdlaed0.f pdlaed1.f pdlaed2.f pdlaed3.f pdlaedz.f 
   pdlamch.f pdlared1d.f pdlasrt.f pdstebz.f pdstedc.f dlapst.f dlasrt2.f 
   dstein2.f
   dlar1va.f dlarrb2.f dlarrd2.f dlarre2.f dlarre2a.f dlarrf2.f
   dlarrv2.f dstegr2.f dstegr2a.f dstegr2b.f
   dbdsdc.f
   dbdsqr.f ddisna.f dlabad.f dlacpy.f dladiv.f dlae2.f  dlaebz.f
   dlaed0.f dlaed1.f dlaed2.f dlaed3.f dlaed4.f dlaed5.f dlaed6.f
   dlaed7.f dlaed8.f dlaed9.f dlaeda.f dlaev2.f dlagtf.f
   dlagts.f dlamrg.f dlanst.f
   dlapy2.f dlapy3.f dlarnv.f
   dlarra.f dlarrb.f dlarrc.f dlarrd.f dlarre.f dlarrf.f dlarrj.f
   dlarrk.f dlarrr.f dlaneg.f
   dlartg.f dlaruv.f dlas2.f  dlascl.f
   dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f
   dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f
   dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f
   dlasr.f  dlasrt.f dlassq.f dlasv2.f dpttrf.f dstebz.f dstedc.f
   dsteqr.f dsterf.f dlaisnan.f disnan.f
   dlartgp.f dlartgs.f
   ../INSTALL/dlamch.f ${DSECOND_SRC})


set (DLASRC  
   pddbsv.f  pddbtrf.f pddbtrs.f pddbtrsv.f 
   pddtsv.f  pddttrf.f pddttrs.f pddttrsv.f 
   pdgbsv.f  pdgbtrf.f pdgbtrs.f 
   pdgebd2.f pdgebrd.f pdgecon.f           pdgeequ.f pdgehd2.f pdgehrd.f 
   pdgelq2.f pdgelqf.f pdgels.f  pdgeql2.f pdgeqlf.f pdgeqpf.f pdgeqr2.f 
   pdgeqrf.f pdgerfs.f pdgerq2.f pdgerqf.f                               
   pdgesv.f  pdgesvd.f pdgesvx.f pdgetf2.f pdgetrf.f pdgetri.f pdgetrs.f 
   pdggqrf.f pdggrqf.f pdlabrd.f 
   pdlacon.f pdlacp2.f pdlacpy.f pdlahrd.f pdlange.f pdlanhs.f pdlansy.f 
   pdlantr.f pdlapiv.f pdlapv2.f pdlaqge.f pdlaqsy.f pdlarf.f  pdlarfb.f 
   pdlarfg.f pdlarft.f pdlase2.f pdlaset.f pdlascl.f pdlassq.f pdlaswp.f 
   pdlatra.f pdlatrd.f pdlatrs.f pdlauu2.f pdlauum.f pdorg2l.f pdorg2r.f 
                       pdorgl2.f pdorglq.f pdorgql.f pdorgqr.f pdorgr2.f 
   pdorgrq.f           pdorm2l.f pdorm2r.f pdormbr.f pdormhr.f pdorml2.f 
   pdormlq.f pdormql.f pdormqr.f pdormr2.f pdormrq.f pdormtr.f pdpocon.f 
   pdpbsv.f  pdpbtrf.f pdpbtrs.f pdpbtrsv.f 
   pdptsv.f  pdpttrf.f pdpttrs.f pdpttrsv.f 
   pdpoequ.f pdporfs.f pdposv.f  pdposvx.f pdpotf2.f pdpotrf.f pdpotri.f 
   pdpotrs.f pdrscl.f  pdstein.f pdsyev.f  pdsyevd.f pdsyevx.f 
   pdsygs2.f pdsygst.f 
   pdsygvx.f pdsyngst.f pdsyntrd.f pdsyttrd.f pdsytd2.f pdsytrd.f pdtrti2.f 
   pdtrtri.f pdtrtrs.f 
   pdlaevswp.f 
   pdlarzb.f pdlarzt.f pdlarz.f pdlatrz.f pdtzrzf.f pdormr3.f pdormrz.f 
   pdlahqr.f pdlaconsb.f pdlacp3.f pdlawil.f 
   pdlasmsub.f pdlared2d.f pdlamr1d.f dlaref.f dlamsh.f dlasorte.f dsteqr2.f 
   ddbtf2.f  ddbtrf.f  ddttrf.f ddttrsv.f dpttrsv.f dtrmvt.f  pdsyevr.f
   bdlaapp.f bdlaexc.f bdtrexc.f dlaqr6.f pdtrord.f
   pdtrsen.f pdgebal.f pdhseqr.f pdlamve.f pdlaqr0.f pdlaqr1.f pdlaqr2.f
   pdlaqr3.f pdlaqr4.f pdlaqr5.f pdrot.f
   dbdsvdx.f dgbbrd.f dgbcon.f dgbequ.f dgbrfs.f dgbsv.f
   dgbsvx.f dgbtf2.f dgbtrf.f dgbtrs.f dgebak.f dgebal.f dgebd2.f
   dgebrd.f dgecon.f dgeequ.f dgees.f  dgeesx.f dgeev.f  dgeevx.f
   dgehd2.f dgehrd.f dgelq2.f dgelqf.f
   dgels.f  dgelsd.f dgelss.f dgelsy.f dgeql2.f dgeqlf.f
   dgeqp3.f dgeqr2.f dgeqr2p.f dgeqrf.f dgeqrfp.f dgerfs.f dgerq2.f dgerqf.f
   dgesc2.f dgesdd.f dgesv.f  dgesvd.f dgesvdx.f dgesvx.f dgetc2.f dgetf2.f
   dgetrf.f dgetrf2.f dgetri.f
   dgetrs.f dggbak.f dggbal.f
   dgges.f  dgges3.f dggesx.f dggev.f  dggev3.f dggevx.f
   dggglm.f dgghrd.f dgghd3.f dgglse.f dggqrf.f
   dggrqf.f dggsvd3.f dggsvp3.f dgtcon.f dgtrfs.f dgtsv.f
   dgtsvx.f dgttrf.f dgttrs.f dgtts2.f dhgeqz.f
   dhsein.f dhseqr.f dlabrd.f dlacon.f dlacn2.f
   dlaein.f dlaexc.f dlag2.f  dlags2.f dlagtm.f dlagv2.f dlahqr.f
   dlahr2.f dlaic1.f dlaln2.f dlals0.f dlalsa.f dlalsd.f
   dlangb.f dlange.f dlangt.f dlanhs.f dlansb.f dlansp.f
   dlansy.f dlantb.f dlantp.f dlantr.f dlanv2.f
   dlapll.f dlapmt.f
   dlaqgb.f dlaqge.f dlaqp2.f dlaqps.f dlaqsb.f dlaqsp.f dlaqsy.f
   dlaqr0.f dlaqr1.f dlaqr2.f dlaqr3.f dlaqr4.f dlaqr5.f
   dlaqtr.f dlar1v.f dlar2v.f iladlr.f iladlc.f
   dlarf.f  dlarfb.f dlarfg.f dlarfgp.f dlarft.f dlarfx.f dlarfy.f
   dlargv.f dlarrv.f dlartv.f
   dlarz.f  dlarzb.f dlarzt.f dlaswp.f dlasy2.f 
   dlasyf.f dlasyf_rook.f dlasyf_rk.f dlasyf_aa.f
   dlatbs.f dlatdf.f dlatps.f dlatrd.f dlatrs.f dlatrz.f dlauu2.f
   dlauum.f dopgtr.f dopmtr.f dorg2l.f dorg2r.f
   dorgbr.f dorghr.f dorgl2.f dorglq.f dorgql.f dorgqr.f dorgr2.f
   dorgrq.f dorgtr.f dorm2l.f dorm2r.f dorm22.f
   dormbr.f dormhr.f dorml2.f dormlq.f dormql.f dormqr.f dormr2.f
   dormr3.f dormrq.f dormrz.f dormtr.f dpbcon.f dpbequ.f dpbrfs.f
   dpbstf.f dpbsv.f  dpbsvx.f
   dpbtf2.f dpbtrf.f dpbtrs.f dpocon.f dpoequ.f dporfs.f dposv.f
   dposvx.f dpotf2.f dpotrf.f dpotrf2.f dpotri.f dpotrs.f dpstrf.f dpstf2.f
   dppcon.f dppequ.f
   dpprfs.f dppsv.f  dppsvx.f dpptrf.f dpptri.f dpptrs.f dptcon.f
   dpteqr.f dptrfs.f dptsv.f  dptsvx.f dpttrs.f dptts2.f drscl.f
   dsbev.f  dsbevd.f dsbevx.f dsbgst.f dsbgv.f  dsbgvd.f dsbgvx.f
   dsbtrd.f dspcon.f dspev.f  dspevd.f dspevx.f dspgst.f
   dspgv.f  dspgvd.f dspgvx.f dsprfs.f dspsv.f  dspsvx.f dsptrd.f
   dsptrf.f dsptri.f dsptrs.f dstegr.f dstein.f dstev.f  dstevd.f dstevr.f
   dstevx.f dsycon.f dsyev.f  dsyevd.f dsyevr.f
   dsyevx.f dsygs2.f dsygst.f dsygv.f  dsygvd.f dsygvx.f dsyrfs.f
   dsysv.f  dsysvx.f
   dsytd2.f dsytf2.f dsytrd.f dsytrf.f dsytri.f dsytrs.f dsytrs2.f
   dsytri2.f dsytri2x.f dsyswapr.f
   dsyconv.f dsyconvf.f dsyconvf_rook.f
   dsytf2_rook.f dsytrf_rook.f dsytrs_rook.f
   dsytri_rook.f dsycon_rook.f dsysv_rook.f
   dsytf2_rk.f dsytrf_rk.f dsytrs_3.f
   dsytri_3.f dsytri_3x.f dsycon_3.f dsysv_rk.f
   dsysv_aa.f dsytrf_aa.f dsytrs_aa.f
   dsysv_aa_2stage.f dsytrf_aa_2stage.f dsytrs_aa_2stage.f
   dtbcon.f
   dtbrfs.f dtbtrs.f dtgevc.f dtgex2.f dtgexc.f dtgsen.f
   dtgsja.f dtgsna.f dtgsy2.f dtgsyl.f dtpcon.f dtprfs.f dtptri.f
   dtptrs.f
   dtrcon.f dtrevc.f dtrevc3.f dtrexc.f dtrrfs.f dtrsen.f dtrsna.f dtrsyl.f
   dtrti2.f dtrtri.f dtrtrs.f dtzrzf.f dstemr.f
   dsgesv.f dsposv.f dlag2s.f slag2d.f dlat2s.f
   dlansf.f dpftrf.f dpftri.f dpftrs.f dsfrk.f dtfsm.f dtftri.f dtfttp.f
   dtfttr.f dtpttf.f dtpttr.f dtrttf.f dtrttp.f
   dgejsv.f dgesvj.f dgsvj0.f dgsvj1.f
   dgeequb.f dsyequb.f dpoequb.f dgbequb.f
   dbbcsd.f dlapmr.f dorbdb.f dorbdb1.f dorbdb2.f dorbdb3.f dorbdb4.f
   dorbdb5.f dorbdb6.f dorcsd.f dorcsd2by1.f
   dgeqrt.f dgeqrt2.f dgeqrt3.f dgemqrt.f
   dtpqrt.f dtpqrt2.f dtpmqrt.f dtprfb.f
   dgelqt.f dgelqt3.f dgemlqt.f
   dgetsls.f dgeqr.f dlatsqr.f dlamtsqr.f dgemqr.f
   dgelq.f dlaswlq.f dlamswlq.f dgemlq.f
   dtplqt.f dtplqt2.f dtpmlqt.f
   dsytrd_2stage.f dsytrd_sy2sb.f dsytrd_sb2st.F dsb2st_kernels.f
   dsyevd_2stage.f dsyev_2stage.f dsyevx_2stage.f dsyevr_2stage.f
   dsbev_2stage.f dsbevx_2stage.f dsbevd_2stage.f dsygv_2stage.f)

set (ZLASRC  
   pzdbsv.f  pzdbtrf.f pzdbtrs.f pzdbtrsv.f 
   pzdtsv.f  pzdttrf.f pzdttrs.f pzdttrsv.f 
   pzgbsv.f  pzgbtrf.f pzgbtrs.f 
   pzgebd2.f pzgebrd.f pzgecon.f           pzgeequ.f pzgehd2.f pzgehrd.f 
   pzgelq2.f pzgelqf.f pzgels.f  pzgeql2.f pzgeqlf.f pzgeqpf.f pzgeqr2.f 
   pzgeqrf.f pzgerfs.f pzgerq2.f pzgerqf.f                               
   pzgesv.f  pzgesvd.f pzgesvx.f pzgetf2.f pzgetrf.f pzgetri.f pzgetrs.f 
   pzggqrf.f 
   pzggrqf.f pzheev.f  pzheevd.f pzheevx.f pzhegs2.f pzhegst.f pzhegvx.f 
   pzhengst.f pzhentrd.f pzhettrd.f pzhetd2.f 
   pzhetrd.f pzlabrd.f pzlacon.f pzlacgv.f pzlacp2.f pzlacpy.f pzlahrd.f 
   pzlahqr.f pzlaconsb.f pzlasmsub.f pzlacp3.f pzlawil.f 
   pzlange.f pzlanhe.f pzlanhs.f pzlansy.f pzlantr.f pzlapiv.f pzlapv2.f 
   pzlaqge.f pzlaqsy.f pzlarf.f  pzlarfb.f pzlarfc.f pzlarfg.f pzlarft.f 
   pzlascl.f pzlase2.f pzlaset.f pzlassq.f pzlaswp.f pzlatra.f pzlatrd.f 
   pzlattrs.f 
   pzlatrs.f pzlauu2.f pzlauum.f pzpocon.f pzpoequ.f pzporfs.f pzposv.f  
   pzpbsv.f  pzpbtrf.f pzpbtrs.f pzpbtrsv.f 
   pzptsv.f  pzpttrf.f pzpttrs.f pzpttrsv.f 
   pzposvx.f pzpotf2.f pzpotrf.f pzpotri.f pzpotrs.f pzdrscl.f pzstein.f 
   pztrevc.f pztrti2.f pztrtri.f pztrtrs.f pzung2l.f pzung2r.f 
             pzungl2.f pzunglq.f pzungql.f pzungqr.f pzungr2.f pzungrq.f 
             pzunm2l.f pzunm2r.f pzunmbr.f pzunmhr.f pzunml2.f pzunmlq.f 
   pzunmql.f pzunmqr.f pzunmr2.f pzunmrq.f pzunmtr.f 
   pzlaevswp.f 
   pzlarzb.f pzlarzt.f pzlarz.f  pzlarzc.f pzlatrz.f pztzrzf.f 
   pzunmr3.f pzunmrz.f pzmax1.f pdzsum1.f pzlamr1d.f 
   zdbtf2.f  zdbtrf.f  zdttrf.f zdttrsv.f  zpttrsv.f zsteqr2.f ztrmvt.f 
   zlamsh.f  zlaref.f  zlanv2.f zlahqr2.f pzheevr.f
   zbdsqr.f zgbbrd.f zgbcon.f zgbequ.f zgbrfs.f zgbsv.f  zgbsvx.f
   zgbtf2.f zgbtrf.f zgbtrs.f zgebak.f zgebal.f zgebd2.f zgebrd.f
   zgecon.f zgeequ.f zgees.f  zgeesx.f zgeev.f  zgeevx.f
   zgehd2.f zgehrd.f zgelq2.f zgelqf.f
   zgels.f  zgelsd.f zgelss.f zgelsy.f zgeql2.f zgeqlf.f zgeqp3.f
   zgeqr2.f zgeqr2p.f zgeqrf.f zgeqrfp.f zgerfs.f zgerq2.f zgerqf.f
   zgesc2.f zgesdd.f zgesv.f  zgesvd.f zgesvdx.f zgesvx.f
   zgesvj.f zgejsv.f zgsvj0.f zgsvj1.f
   zgetc2.f zgetf2.f zgetrf.f zgetrf2.f
   zgetri.f zgetrs.f
   zggbak.f zggbal.f
   zgges.f  zgges3.f zggesx.f zggev.f  zggev3.f zggevx.f
   zggglm.f zgghrd.f zgghd3.f zgglse.f zggqrf.f zggrqf.f
   zggsvd3.f zggsvp3.f
   zgtcon.f zgtrfs.f zgtsv.f  zgtsvx.f zgttrf.f zgttrs.f zgtts2.f zhbev.f
   zhbevd.f zhbevx.f zhbgst.f zhbgv.f  zhbgvd.f zhbgvx.f zhbtrd.f
   zhecon.f zheev.f  zheevd.f zheevr.f zheevx.f zhegs2.f zhegst.f
   zhegv.f  zhegvd.f zhegvx.f zherfs.f zhesv.f  zhesvx.f zhetd2.f
   zhetf2.f zhetrd.f
   zhetrf.f zhetri.f zhetri2.f zhetri2x.f zheswapr.f
   zhetrs.f zhetrs2.f
   zhetf2_rook.f zhetrf_rook.f zhetri_rook.f
   zhetrs_rook.f zhecon_rook.f zhesv_rook.f
   zhetf2_rk.f zhetrf_rk.f zhetri_3.f zhetri_3x.f
   zhetrs_3.f zhecon_3.f zhesv_rk.f
   zhesv_aa.f zhetrf_aa.f zhetrs_aa.f
   zhesv_aa_2stage.f zhetrf_aa_2stage.f zhetrs_aa_2stage.f
   zhgeqz.f zhpcon.f zhpev.f  zhpevd.f
   zhpevx.f zhpgst.f zhpgv.f  zhpgvd.f zhpgvx.f zhprfs.f zhpsv.f
   zhpsvx.f
   zhptrd.f zhptrf.f zhptri.f zhptrs.f zhsein.f zhseqr.f zlabrd.f
   zlacgv.f zlacon.f zlacn2.f zlacp2.f zlacpy.f zlacrm.f zlacrt.f zladiv.f
   zlaed0.f zlaed7.f zlaed8.f
   zlaein.f zlaesy.f zlaev2.f zlags2.f zlagtm.f
   zlahef.f zlahef_rook.f zlahef_rk.f zlahef_aa.f zlahqr.f
   zlahr2.f zlaic1.f zlals0.f zlalsa.f zlalsd.f zlangb.f zlange.f
   zlangt.f zlanhb.f
   zlanhe.f
   zlanhp.f zlanhs.f zlanht.f zlansb.f zlansp.f zlansy.f zlantb.f
   zlantp.f zlantr.f zlapll.f zlapmt.f zlaqgb.f zlaqge.f
   zlaqhb.f zlaqhe.f zlaqhp.f zlaqp2.f zlaqps.f zlaqsb.f
   zlaqr0.f zlaqr1.f zlaqr2.f zlaqr3.f zlaqr4.f zlaqr5.f
   zlaqsp.f zlaqsy.f zlar1v.f zlar2v.f ilazlr.f ilazlc.f
   zlarcm.f zlarf.f  zlarfb.f
   zlarfg.f zlarfgp.f zlarft.f
   zlarfx.f zlarfy.f zlargv.f zlarnv.f zlarrv.f zlartg.f zlartv.f
   zlarz.f  zlarzb.f zlarzt.f zlascl.f zlaset.f zlasr.f
   zlassq.f zlaswp.f zlasyf.f zlasyf_rook.f zlasyf_rk.f zlasyf_aa.f
   zlatbs.f zlatdf.f zlatps.f zlatrd.f zlatrs.f zlatrz.f zlauu2.f
   zlauum.f zpbcon.f zpbequ.f zpbrfs.f zpbstf.f zpbsv.f
   zpbsvx.f zpbtf2.f zpbtrf.f zpbtrs.f zpocon.f zpoequ.f zporfs.f
   zposv.f  zposvx.f zpotf2.f zpotrf.f zpotrf2.f zpotri.f zpotrs.f zpstrf.f zpstf2.f
   zppcon.f zppequ.f zpprfs.f zppsv.f  zppsvx.f zpptrf.f zpptri.f zpptrs.f
   zptcon.f zpteqr.f zptrfs.f zptsv.f  zptsvx.f zpttrf.f zpttrs.f zptts2.f
   zrot.f   zspcon.f zspmv.f  zspr.f   zsprfs.f zspsv.f
   zspsvx.f zsptrf.f zsptri.f zsptrs.f zdrscl.f zstedc.f
   zstegr.f zstein.f zsteqr.f zsycon.f zsymv.f
   zsyr.f   zsyrfs.f zsysv.f  zsysvx.f zsytf2.f zsytrf.f zsytri.f
   zsytri2.f zsytri2x.f zsyswapr.f
   zsytrs.f zsytrs2.f
   zsyconv.f zsyconvf.f zsyconvf_rook.f
   zsytf2_rook.f zsytrf_rook.f zsytrs_rook.f zsytrs_aa.f
   zsytri_rook.f zsycon_rook.f zsysv_rook.f
   zsytf2_rk.f zsytrf_rk.f zsytrf_aa.f zsytrs_3.f
   zsysv_aa_2stage.f zsytrf_aa_2stage.f zsytrs_aa_2stage.f
   zsytri_3.f zsytri_3x.f zsycon_3.f zsysv_rk.f zsysv_aa.f
   ztbcon.f ztbrfs.f ztbtrs.f ztgevc.f ztgex2.f
   ztgexc.f ztgsen.f ztgsja.f ztgsna.f ztgsy2.f ztgsyl.f ztpcon.f
   ztprfs.f ztptri.f
   ztptrs.f ztrcon.f ztrevc.f ztrevc3.f ztrexc.f ztrrfs.f ztrsen.f ztrsna.f
   ztrsyl.f ztrti2.f ztrtri.f ztrtrs.f ztzrzf.f zung2l.f
   zung2r.f zungbr.f zunghr.f zungl2.f zunglq.f zungql.f zungqr.f zungr2.f
   zungrq.f zungtr.f zunm2l.f zunm2r.f zunmbr.f zunmhr.f zunml2.f zunm22.f
   zunmlq.f zunmql.f zunmqr.f zunmr2.f zunmr3.f zunmrq.f zunmrz.f
   zunmtr.f zupgtr.f
   zupmtr.f izmax1.f dzsum1.f zstemr.f
   zcgesv.f zcposv.f zlag2c.f clag2z.f zlat2c.f
   zhfrk.f ztfttp.f zlanhf.f zpftrf.f zpftri.f zpftrs.f ztfsm.f ztftri.f
   ztfttr.f ztpttf.f ztpttr.f ztrttf.f ztrttp.f
   zgeequb.f zgbequb.f zsyequb.f zpoequb.f zheequb.f
   zbbcsd.f zlapmr.f zunbdb.f zunbdb1.f zunbdb2.f zunbdb3.f zunbdb4.f
   zunbdb5.f zunbdb6.f zuncsd.f zuncsd2by1.f
   zgeqrt.f zgeqrt2.f zgeqrt3.f zgemqrt.f
   ztpqrt.f ztpqrt2.f ztpmqrt.f ztprfb.f
   ztplqt.f ztplqt2.f ztpmlqt.f
   zgelqt.f zgelqt3.f zgemlqt.f
   zgetsls.f zgeqr.f zlatsqr.f zlamtsqr.f zgemqr.f
   zgelq.f zlaswlq.f zlamswlq.f zgemlq.f
   zhetrd_2stage.f zhetrd_he2hb.f zhetrd_hb2st.F zhb2st_kernels.f
   zheevd_2stage.f zheev_2stage.f zheevx_2stage.f zheevr_2stage.f
   zhbev_2stage.f zhbevx_2stage.f zhbevd_2stage.f zhegv_2stage.f)
   
set(src 
      ${ALLAUX} ${DZLAUX}
      ${DLASRC} ${ZLASRC})

set(src-C 
      ${ALLAUX-C} )

################################################################3



