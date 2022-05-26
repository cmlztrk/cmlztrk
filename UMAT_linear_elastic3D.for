
!-----------------Start of Base code(UMAT_ISOTROP_3D), DO NOT CHANGE ANYTHING!!!---------

	  SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
!

      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
!--------------------------------------------------------------------------------
!---------It is at the END OF BASE CODE------------------------------------------

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

!--------------------Start of USER code------------------------------------------
!--------------------------------------------------------------------------------
!              PROPS(NPROPS)
!              User-specified array of material constants associated with this user material.
! 
!              NPORPS
!              NPROPS
!              User-defined number of material constants associated with this user material.
! 
!              PROPS(1) = Elasticity modulus
!              PROPS(2) = Poisson`s Ratio
!              G = Shear Modulus
		E=PROPS(1)
		v=PROPS(2)
		G=E/(2.D0*(1.D0+v))
!--------------------------------------------------------------------------------
!              DDSDDE(NTENS,NTENS)
!              Jacobian matrix or strain component array (NDI + NSHR)

!              NDI
!              Number of direct stress components at this point
!--------------------------------------------------------------------------------
!              NSHR
!              Number of engineering shear stress components at this point.

		DO i=1, NDI
			DO j=1, NDI
				DDSDDE(j, i)=(E*v)/((1+v)*(1-2*v))
			END DO
			DDSDDE(i, i)=(E*(1-v))/((1+v)*(1-2*v))
		END DO
	
	
		DO i=NDI+1, NTENS
			DDSDDE(i ,i)=G
		END DO
	
!---------------Calculation of stresses------------------------------------------
!--------------------------------------------------------------------------------
!               STRESS(NTENS)
!               This array is passed in as the stress tensor at the beginning of the oncrement and must be updated
!               in this routine to be the stress tensor at the end of the increment.
!
!               DSTRAN(NTENS)
!               Array of strain increments.

		DO i=1, NTENS
		DO j=1, NTENS
			STRESS(j)=STRESS(j)+DDSDDE(j, i)*DSTRAN(i)
		END DO
		END DO
!	
!-------------End of USER code--------------------------------------------------
!
		RETURN
		END