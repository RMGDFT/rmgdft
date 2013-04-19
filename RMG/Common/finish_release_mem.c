/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"

void finish_release_mem (STATE * states)
{


    int ion, isp;
    ION *iptr;
    SPECIES *sp;




    if (pct.weight)
        my_free (pct.weight);

    /*************** Memory for non-local projectors **************************/
    for (ion = 0; ion < ct.num_ions; ion++)
    {


#if FDIFF_BETA

        if (pct.weight_derx[ion])
            my_free (pct.weight_derx[ion]);
        if (pct.weight_dery[ion])
            my_free (pct.weight_dery[ion]);
        if (pct.weight_derz[ion])
            my_free (pct.weight_derz[ion]);

#endif

        if (pct.nlindex[ion])
            my_free (pct.nlindex[ion]);

        if (pct.Qindex[ion])
            my_free (pct.Qindex[ion]);

        if (pct.idxflag[ion])
            my_free (pct.idxflag[ion]);

        if (pct.Qdvec[ion])
            my_free (pct.Qdvec[ion]);



        if (pct.phaseptr[ion])
            my_free (pct.phaseptr[ion]);

        if (pct.augfunc[ion])
            my_free (pct.augfunc[ion]);

        if (pct.dnmI[ion])
            my_free (pct.dnmI[ion]);

        if (pct.qqq[ion])
            my_free (pct.qqq[ion]);

    }





#if FDIFF_BETA
    if (pct.weight_derx)
        my_free (pct.weight_derx);

    if (pct.weight_dery)
        my_free (pct.weight_dery);

    if (pct.weight_derz)
        my_free (pct.weight_derz);
#endif

    if (pct.nlindex)
        my_free (pct.nlindex);

    if (pct.Qindex)
        my_free (pct.Qindex);

    if (pct.idxflag)
        my_free (pct.idxflag);

    if (pct.Qdvec)
        my_free (pct.Qdvec);

    if (pct.idxptrlen)
        my_free (pct.idxptrlen);

    if (pct.Qidxptrlen)
        my_free (pct.Qidxptrlen);

    if (pct.lptrlen)
        my_free (pct.lptrlen);

    if (pct.phaseptr)
        my_free (pct.phaseptr);

    if (pct.augfunc)
        my_free (pct.augfunc);

    if (pct.dnmI)
        my_free (pct.dnmI);

    if (pct.qqq)
        my_free (pct.qqq);
    /*************** End for non-local projectors **************************/



    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {


        /* Get species type */
        sp = &ct.sp[isp];

        /*Memory from init_derweight */
        if (sp->forward_derbeta_x)
            my_free (sp->forward_derbeta_x);


        /*Memory from init_weight */
        if (sp->forward_beta)
            my_free (sp->forward_beta);


        /*Memory from read_pseudo */
        if (sp->qfcoef)
            my_free (sp->qfcoef);

        if (sp->qnm)
            my_free (sp->qnm);

        if (sp->qnmlig)
            my_free (sp->qnmlig);

        if (sp->drqnmlig)
            my_free (sp->drqnmlig);


    }


    /*Fourier transform phase */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];

        if (iptr->fftw_phase_sin)
            my_free (iptr->fftw_phase_sin);

        if (iptr->fftw_phase_cos)
            my_free (iptr->fftw_phase_cos);

    }


    /*Memory for wavefunctions */
    if (states[0].psiR)
        my_free (states[0].psiR);



    if (ct.vh_ext)
        my_free (ct.vh_ext);


    if (pct.newsintR_local) my_free (pct.newsintR_local);
    if (pct.oldsintR_local) my_free (pct.oldsintR_local);

#if !GAMMA_PT
    if (pct.newsintI_local) my_free (pct.newsintI_local);
    if (pct.oldsintI_local) my_free (pct.oldsintI_local);
#endif


    if (states)
        my_free (states);



    /*Release memory for species */
    if (ct.sp)
        my_free (ct.sp);


    /*Release memory for ions */
    if (ct.ions)
        my_free (ct.ions);


    /*Release memory for k-point structure */
    if (ct.kp)
        my_free (ct.kp);




}

/*EOF*/
