
void mkS_canonic_interval(cst_pI  	degree,
			   cst_pI 	n,
			   pR 	r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   pR 	rwork,
			   cst_pI 	rwork_n,
			   pI 	err_)
{
  err_[0] = (I)0;
  const I k = degree[0];
  switch(k)
    {
    case 0:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	    } }
	break;
      }
    case 1:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	      r[roff_[0]*i+1] = p[poff_[0]*i];

	    } }
	break;
      }
    case 2:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	      r[roff_[0]*i+1] = p[poff_[0]*i];
	      r[roff_[0]*i+2] = p[poff_[0]*i]*p[poff_[0]*i];
	    } }
	break;
      }
    case 3:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	      r[roff_[0]*i+1] = p[poff_[0]*i];
	      r[roff_[0]*i+2] = p[poff_[0]*i]*p[poff_[0]*i];
	      r[roff_[0]*i+3] = p[poff_[0]*i]*r[roff_[0]*i+2];
	    } }
	break;
      }
    default:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	      r[roff_[0]*i+1] = p[poff_[0]*i];
	      r[roff_[0]*i+2] = p[poff_[0]*i]*p[poff_[0]*i];
	      r[roff_[0]*i+3] = p[poff_[0]*i]*r[roff_[0]*i+2];
	      { I j;
		for (j=4;j<=k;++j)
		  {
		    r[roff_[0]*i+j] = p[poff_[0]*i]*r[roff_[0]*i+(j-1)];
		  } }
	    } }
	break;
      }
    }
}



void mkS_dx_canonic_interval(cst_pI  	degree,
			     cst_pI 	n,
			     pR 	r,
			     cst_pI 	roff_,
			     cst_pR 	p,
			     cst_pI 	poff_,
			     pR 	rwork,
			     cst_pI 	rwork_n,
			     pI 	err_)
{
  err_[0] = (I)0;
  const I k = degree[0];
  switch(k)
    {
    case 0:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	    } }
	break;
      }
    case 1:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	      r[roff_[0]*i+1] = ((R)1.0);

	    } }
	break;
      }
    case 2:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	      r[roff_[0]*i+1] = ((R)1.0);
	      r[roff_[0]*i+2] = ((R)2.0)*p[poff_[0]*i];
	    } }
	break;
      }
    case 3:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	      r[roff_[0]*i+1] = ((R)1.0);
	      r[roff_[0]*i+2] = ((R)2.0)*p[poff_[0]*i];
	      r[roff_[0]*i+3] = ((R)3.0)*p[poff_[0]*i]*p[poff_[0]*i];
	    } }
	break;
      }
    default:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	      r[roff_[0]*i+1] = ((R)1.0);
	      r[roff_[0]*i+2] = ((R)2.0)*p[poff_[0]*i];
	      r[roff_[0]*i+3] = ((R)3.0)*p[poff_[0]*i]*p[poff_[0]*i];
	      { I j;
		for (j=4;j<=k;++j)
		  {
		    r[roff_[0]*i+j] = p[poff_[0]*i]*( r[roff_[0]*i+(j-1)]/((R)(j-1)) ) * ((R)j);
		  } }
	    } }
	break;
      }
    }
}
