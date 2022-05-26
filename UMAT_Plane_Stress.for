!-----------------This code was created thanks to Ronald Wagner - You can check this link for more information.
!                       (https://www.youtube.com/watch?v=zDeb9WGsC-8&ab_channel=Dr.-Ing.RonaldWagner)
!-----------------Start of Base code(UMAT_ISOTROP_3D), DO NOT CHANGE ANYTHING!!!-------------------------------
!
!
!
!-----------------Start of Base code(UMAT_ISOTROP_3D), DO NOT CHANGE ANYTHING!!!-------------------------------
	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

!--------------It is the at the end of BASE code UMAT--------------------------------------------------------
!
		PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
!
!--------------It is the part of USER code--------------------------------------------------------------------
!              PROPS(NPROPS)
!              User-specified array of material constants associated with this user material.
! 
!              NPORPS
!              NPROPS
!              User-defined number of material constants associated with this user material.
! 
!              PROPS(1) = Elasticity modulus
!              PROPS(2) = Poisson`s Ratio
!
			E=PROPS(1)
			v=PROPS(2)
			
				DO K1 = 1, NTENS
					DO K2 = 1, NTENS
					DDSDDE(K1, K2)=ZERO
					END DO
				END DO
!
!              DDSDDE(NTENS,NTENS)
!              Jacobian matrix or strain component array (NDI + NSHR)
!
!              NDI
!              Number of direct stress components at this point
!--------------------------------------------------------------------------------
!              NSHR
!              Number of engineering shear stress components at this point.
!
!				Stiffness Matrix for Plane Stress
!
				DDSDDE(1,1) = E/(ONE-v*v)
				DDSDDE(2,2) = E/(ONE-v*v)
				DDSDDE(3,3) = E*(ONE-v)/(TWO*(ONE-v*v))
					
				DDSDDE(1,2) = E*v/(ONE-v*v)
				DDSDDE(2,1) = E*v/(ONE-v*v)
!
!---------------Calculation of stresses------------------------------------------
!--------------------------------------------------------------------------------
!               STRESS(NTENS)
!               This array is passed in as the stress tensor at the beginning of the oncrement and must be updated
!               in this routine to be the stress tensor at the end of the increment.
!
!               DSTRAN(NTENS)
!               Array of strain increments.
!

				DO K1=1, NTENS
					DO K2=1, NTENS
					STRESS(K2)=STRESS(K2)+DDSDDE(K1, K2)*DSTRAN(K1)
					END DO
				END DO
				
!--------------------------------------------------------------------------------
!							End of USER code
!--------------------------------------------------------------------------------
	
      RETURN
      END